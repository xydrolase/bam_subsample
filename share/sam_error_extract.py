#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
sam_error_extract.py

Author: Xin Yin <xinyin at iastate dot edu>

This is a postprocessing script for short-read alignment results in SAM format.
If you are using an aligner that can output multiple formats, it is important
to specify the output format as SAM.

Please refer to the SAM Format Specification for more details on 
manipulation of SAM format datafiles.
    http://samtools.sourceforge.net/SAM1.pdf 


After aligning a sequencing dataset to its reference, which yields a SAM file
ready to process, one can use this tool to extract "ground truth" errors 
(all mismatches against the reference genome) for error correction 
benchmarking. 

This tool does not extract ALL errors. Instead, it only tries to extract errors
from a white list of read IDs provided, in the arguments, by the user.
However, not all reads in the white list are eligible for error extraction.
Specifically, reads will be disqualified for error extraction in any of the
following cases:

    - Ambiguous mapping (mapping to multiple genome locations)
    - Containing soft/hard clipping (so that errors are censored)
    - Containing insertion/deletions (so that error position will be shifted)

!!IMPORTANT!!
-------------

All reads that are skipped in the process will be tabulated in a PREFIX.skipped
file. Later, this file can be provided to the benchmarking script so that all
errors originated from these reads will be discarded for benchmarking, avoiding
counting error corrections on these reads as false positives. 
"""

from argparse import ArgumentParser
from itertools import izip, cycle
from math import floor

import sys
import re

complements = {'A':'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
# map cigar operators to how they would affect length of the mapped
# segment.
cigar_op_map = {'M': 1, 'I': 0, 'D': -1, 'S': 1, 'H': 0}

# column definitions
COL_QNAME, COL_FLAG, COL_RNAME, COL_POS, COL_MAPQ, COL_CIGAR, \
        COL_RNEXT, COL_PNEXT, COL_TLEN, COL_SEQ, COL_QUAL = range(11)

# mapping flags
FLAG_UNMAPPED = 0x4
FLAG_REVCOMP = 0x10
FLAG_1ST_READ = 0x40
FLAG_2ND_READ = 0x80
FLAG_SECONDARY = 0x100

def compute_segment_length(c_op, c_num):
    return sum([x*y for x, y in zip(c_num, 
        [cigar_op_map.get(op, 0) for op in c_op])])

def compute_base_offset(c_op, c_num):
    if c_op[0] == 'S':
        return c_num[0]

    return 0

def main():
    parser = ArgumentParser(
        description="Genome sequence guided FASTQ file subsampler.")
    parser.add_argument("-Q", dest="output_qscore", action="store_true",
            default=False, 
            help="Enable output of quality scores in error file.")
    parser.add_argument("--prefix", dest="prefix", required=True,
            help="Prefix for output files.")
    parser.add_argument("-C", "--skip-clipping", dest="skip_clipping",
            action="store_true", default=False,
            help="Skip all the reads with soft/hard clippings.")
    parser.add_argument("-l", "--length", metavar="L", dest="length",
            type=int, help="Cap the maximum read length to L.")
    parser.add_argument("--discard-N", dest="discard_N", 
            action="store_true", default=False,
            help="Skip all the reads with N base(s).")
    parser.add_argument("--expand-read-id", dest="expand_read_id",
            action="store_true", default=False,
            help="For SAM file generated by Novoalign, expand read id to "
            "include matepair information.")
    
    parser.add_argument("read_id_list", type=file, 
            help="Sequencing reads in FASTQ format.")
    parser.add_argument("sam_file", type=file,
            help="Alignment results in SAM format.")

    args = parser.parse_args()

    fd_error = open('{0}.error'.format(args.prefix), 'w+')
    fd_skip = open('{0}.skipped'.format(args.prefix), 'w+')

    regex_base = re.compile(r'([ACGT])')
    regex_cigar = re.compile(r'[MIDNSH]')
    regex_clip = re.compile(r'\d+[SH]')
    set_bases = set(['A', 'C', 'G', 'T'])

    selected_reads = set([line.strip() for line in args.read_id_list])

    for line in args.sam_file:
        if not selected_reads:
            break

        # skip headers
        if line.startswith('@'):
            continue

        segments = line.strip().split('\t')

        align_flag = int(segments[COL_FLAG])
        map_qual = int(segments[COL_MAPQ])
        map_pos = int(segments[COL_POS])
        cigar = segments[COL_CIGAR]

        # process errors
        read_id = segments[COL_QNAME]

        if args.expand_read_id:
            if align_flag & FLAG_1ST_READ:
                read_id = '%s/1' % read_id
            elif align_flag & FLAG_2ND_READ:
                read_id = '%s/2' % read_id

        if read_id not in selected_reads:
            continue

        selected_reads.remove(read_id)

        # skip all ambiguously mapped/unmapped reads.
        if align_flag & FLAG_UNMAPPED or map_qual == 0 or cigar == '*':
            fd_skip.write(read_id + '\n')
            continue

        # process cigar
        cigar_ops = re.findall(regex_cigar, cigar)
        cigar_nums = [int(x) for x in re.split(regex_cigar, cigar) if x]

        # skip read with soft/hard clippings if user opted so.
        if args.skip_clipping and regex_clip.search(cigar): 
            fd_skip.write(read_id + '\n')
            continue

        # the true segment length in reference genome
        segment_len = compute_segment_length(cigar_ops, cigar_nums)

        read_seq = segments[COL_SEQ]
        read_qscore = segments[COL_QUAL]

        md_tag = [seg.split(':')[2] for seg in segments if seg.startswith('MD')][0]
        if md_tag.find('^') != -1:
            fd_skip.write(read_id + '\n')
            continue

        md_record = re.split(regex_base, md_tag)

        # skip all the reads with insertions/deletions and hard clippings.
        if segment_len != len(read_seq):
            fd_skip.write(read_id + '\n')
            continue

        if args.discard_N and read_seq.find('N') != -1:
            fd_skip.write(read_id + '\n')
            continue

        # find the the base offset to compute the actual error loci
        # for reads with soft clipping.
        read_pos = compute_base_offset(cigar_ops, cigar_nums)

        for md in md_record: 
            if md.isdigit():
                read_pos += int(md)
            elif md in set_bases:
                if align_flag & FLAG_REVCOMP:
                    fd_error.write('{0} {1} {2} {3}{4}\n'.format(
                        read_id, segment_len-read_pos, # 1-based coordinate
                        complements[md], # compute the complementary base
                        complements[read_seq[read_pos]],
                        ' '.join(('', read_qscore[read_pos])) if \
                                args.output_qscore else ''
                        ))
                else:
                    fd_error.write('{0} {1} {2} {3}{4}\n'.format(
                        read_id, read_pos+1,
                        md, read_seq[read_pos],
                        ' '.join(('', read_qscore[read_pos])) if \
                                args.output_qscore else ''
                        ))

                # each mismatch record also accounts for 1 base 
                read_pos += 1


    fd_error.close()
    fd_skip.close()

if __name__ == "__main__":
    main()