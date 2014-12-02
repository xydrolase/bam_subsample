/*  bam_subsample.cpp 
 *  Generating FASTQ and ERROR files from alignment results in BAM format.
 *
 *
 *    Author: Xin Yin <xinyin at iastate dot edu>
 *
 */

#include <cstdlib>
#include <fstream>
#include <boost/program_options.hpp>

extern "C" {
#include "htslib/sam.h" // alignment
#include "htslib/bgzf.h"  // BAM
}

#include "bam_subsample.h"

namespace progopt = boost::program_options;

static inline void format_error(std::ostream &ss, uint32_t is_rev, std::string rid, 
		int rlen, int epos, char ref, char err);
void tally_errors(char *MD, std::vector<std::string> &errs, 
		std::string &rid, const uint8_t *rd_seq, uint32_t flag, int qlen);

// convert 4-bit repr. of a nucleotide to the "true" literal nucleotide ACTG.
static inline char bam_4bit2char(uint8_t b4) 
{
	// use GCC/GNU builtin to scan the number of trailing zeros.
	// i.e. 1 -> 0, 2 -> 1, 4 -> 2, 8 -> 3
	int idx = __builtin_ctz(b4); 
	uint8_t mask = (b4 == 0xf) - 1; // 0xff if ACTG, 0 if N
	// is it the "N" base?
	int nidx = (idx & mask) + (4 & ~mask);
	return sam_nucleotides[nidx];
}

// convert 4-bit repr. of a nucleotide to the "true" literal nucleotide ACTG.
static inline char bam_4bit2comp(uint8_t b4) 
{
	// use GCC/GNU builtin to scan the number of trailing zeros.
	// i.e. 1 -> 0, 2 -> 1, 4 -> 2, 8 -> 3
	// then, evaluate the complementary base.
	// Original:      ACGT 0, 1, 2, 3
	// Complementary: TGCA 0, 1, 2, 3
	int idx = 3 - __builtin_ctz(b4); 
	uint8_t mask = (b4 == 0xf) - 1; // 0xff if ACTG, 0 if N
	// is it the "N" base?
	int nidx = (idx & mask) + (4 & ~mask);
	return sam_nucleotides[nidx];
}

void format_error(std::stringstream &ss, uint32_t is_rev, std::string rid, 
		int qlen, int epos, char ref, char err)
{
	// Format for output:
	// READ_ID POS_IN_READ TRUE_BASE ERROR
	
	if (is_rev) {
		epos = qlen - epos;
		ref = pmr_complements[base_to_bits(ref)];
		err = pmr_complements[base_to_bits(err)];
	}
	else {
		++epos; // convert from 0-base to 1-based coordinate
	}

	ss << rid << " " << epos << " " << ref << " " << err << std::endl;
}

// process MD aux tag
void tally_errors(char *MD, std::vector<std::string> &errs, 
		std::string &rid, const uint8_t *rd_seq, uint32_t flag, int qlen)
{
	uint32_t epos = 0;
	uint32_t radix = 0;

	for (char *p = MD; *p != '\0'; ++p) {
		char pc = *p;
		if (isalpha(pc)) {
			epos += radix;
			// observed nucleotide on the read:
			// NOTE that bam_seqi returns a 4-bit representation of
			// nucleotide: A = 1, C = 2, G = 4, T = 8, N = 15.
			// We need to convert the 4-bit repr. back to char
			char obs_epos_base = bam_4bit2char(bam_seqi(rd_seq, epos));

#ifdef DEBUG
			std::cerr << "Base @ " << epos << " == " << (int) (bam_seqi(rd_seq, epos)) <<
				"(4-bit repr.) -> " << obs_epos_base << std::endl;
#endif
			std::stringstream serr;
			format_error(serr, flag & BAM_FREVERSE, rid, qlen, epos, 
					pc, obs_epos_base);
			errs.push_back(serr.str());
			++epos;

			radix = 0; // reset radix
		}
		else if (isdigit(pc)) {
			radix = radix * 10 + (pc - '0');
		}
		else {
			// likely a ^ char indicating insertion errors
			// TODO: add support for handling insertion errors.
			std::cerr << "Invalid character: " << pc << std::endl;
			exit(1);
		}
	}
}

// A single-threaded BAM file processor, that generates a FASTQ file 
// corresponding to the reads mapped, and an ERROR file reporting the
// sequencing errors (but somehow could be true genomic variants as well).
// Filter options are defined in the `opt` struct.

void bam_filter(options &opt)
{
	BGZF *bam_fp = bgzf_open(opt.bam_file.c_str(), "r");
	bam_hdr_t *hdr = bam_hdr_read(bam_fp);

	std::vector<int> tid_fltr;
	
	// check --reference filter
	for (int i = 0; i < hdr->n_targets; ++i) {
		for (auto it = opt.references.begin(); it != opt.references.end(); ++it) {
			if (it->compare(hdr->target_name[i]) == 0) {
				tid_fltr.push_back(i);
				std::cerr << "Adding filter: " << hdr->target_name[i] << std::endl;
			}
		}
	}

	bam1_t *aln = bam_init1();

	// open output files
	std::ofstream of_fastq;
	std::ofstream of_error;
	std::string _fn = opt.output_prefix + ".fastq";
	of_fastq.open(_fn.c_str(), std::ofstream::out | std::ofstream::trunc);
	_fn = opt.output_prefix + ".error";
	of_error.open(_fn.c_str(), std::ofstream::out | std::ofstream::trunc);

	std::stringstream mtpair_ss;
	std::vector<std::string> mtpair_errors;

	int nbytes = 0;
	int nread = 0;
	bool pmapped = false;
	while ((nbytes = bam_read1(bam_fp, aln)) >= 0) {
		if (nread++ % 1000000 == 0) {
			std::cerr << (nread-1) << " reads processed." << std::endl;
		}

		uint32_t flag = aln->core.flag;

		if (flag & BAM_FREAD1) {
			// reset matepair reads and errors
			mtpair_ss.str(std::string());
			mtpair_errors.clear();
		}

		if (opt.paired_only && (flag & BAM_FREAD2) && !pmapped) {
			pmapped = false;
			continue;
		}
		// reset the flag indicating if the previous read is mapped or not
		pmapped = false;

		// check if unmapped, if so, discard immediately.
		// also discards alignment that is not primary.
		if (flag & BAM_FUNMAP || flag & BAM_FSECONDARY) {
#ifdef DEBUG
			std::cout << bam_get_qname(aln) << std::endl;
			getchar();
#endif
			continue;
		}

		if (!tid_fltr.empty()) {
			int tid = aln->core.tid;
			bool found = false;
			for (auto it = tid_fltr.begin(); it != tid_fltr.end(); ++it) {
				found |= (*it == tid);
			}

			if (!found) continue;
		}

		// expand read identifier according to the flag BAM_FREAD1 / BAM_FREAD2
		
		std::stringstream erid;
		erid << bam_get_qname(aln);
		if (opt.expand_read_id) {
			erid << READ_EXPAND_SEP << ((flag & BAM_FREAD1) ? "1" : "2");
		}

		// retrieve cigar, used to check 
		uint32_t *cigar = bam_get_cigar(aln);
		int n_cigar = aln->core.n_cigar;

		// qlen = Query sequence LENgth
		int qlen = aln->core.l_qseq;

		// length of mapped reference
		int rlen = bam_cigar2rlen(n_cigar, cigar);

		bool has_clip = false;	
		bool has_indel = false;
		for (int i = 0; i < n_cigar; ++i) {
			register uint32_t op = bam_cigar_op(cigar[i]);
			has_clip |= (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP);
			has_indel |= (op == BAM_CINS || op == BAM_CDEL);
		}

		// check clip
		if (opt.skip_clip && has_clip) continue;

		// check indel
		if (opt.skip_indel && has_indel) {
			continue;
		}

		// pointer to the 4-bit repr. of the read sequence
		const uint8_t *rd_seq = bam_get_seq(aln);
		// quality score (literal, without the 33 offset). 
		const char *rd_qual = (const char *) bam_get_qual(aln); 

		char rseq_literal[qlen + 1];
		char rseq_qual[qlen + 1];

		rseq_literal[qlen] = '\0';
		rseq_qual[qlen] = '\0';

		uint8_t contains_n = 0;
		if (flag & BAM_FREVERSE) {
			// for reads that map to the reverse complemnt of the reference,
			// we need to flip the sequence order accordingly (and of course
			// evaluate the complementary sequence).
			for (int i = 0; i < qlen; ++i) {
				int ri = qlen - i - 1; // 0-based index
				uint8_t b = bam_seqi(rd_seq, i);
				rseq_literal[ri] = bam_4bit2comp(b);
				contains_n |= (b == 0xf);

				// convert quality score to ASCII representation. 
				rseq_qual[ri] = rd_qual[i] + 33;
			}
		}
		else {
			for (int i = 0; i < qlen; ++i) {
				uint8_t b = bam_seqi(rd_seq, i);
				rseq_literal[i] = bam_4bit2char(b);
				contains_n |= (b == 0xf);

				// convert quality score to ASCII representation. 
				rseq_qual[i] = rd_qual[i] + 33;
			}
		}

		// check if the sequence contain ambiguous bases
		if (opt.skip_nbase && contains_n) {
			continue;
		}

		// tally mismatches
		uint8_t *aux = bam_aux_get(aln, aux_md);
		char *MD = bam_aux2Z(aux);

		std::string rd_identif = erid.str();
		// FIXME: currently only handles substitution errors,
		// will crash if insertion error is encountered.
		// therefore, always run with -I (skip indel) option. 
		tally_errors(MD, mtpair_errors, rd_identif, rd_seq, flag, qlen);

		// FASTQ record
		mtpair_ss << "@" << rd_identif << std::endl 
			<< rseq_literal << std::endl
			<< "+" << rd_identif << std::endl
			<< rseq_qual << std::endl;

		pmapped = true;

		// require paired?
		if (opt.paired_only && (flag & BAM_FREAD1)) continue;

		// FIXME: output has some bug if -P option is not specified.
		of_fastq << mtpair_ss.str();
		for (auto it = mtpair_errors.begin(); it != mtpair_errors.end(); ++it) {
			of_error << *it;
		}
	}

	of_fastq.close();
	of_error.close();
}

void parse_options(int argc, char *argv[], options &opt)
{
	progopt::options_description desc("BAM filter options");
	desc.add_options()
	("help,h", "print usage information")
	("output,o", progopt::value(&opt.output_prefix), "prefix for output files")
	("reference", progopt::value<std::vector<std::string>>()->composing(), 
	 "reference (template) name for filtering")
	("skip_indel,I", progopt::bool_switch(&opt.skip_indel), 
	 "skip reads containing indel errors")
	("skip_nbase,N", progopt::bool_switch(&opt.skip_nbase),
	 "skip reads containing ambiguous (N) bases")
	("skip_clip,C", progopt::bool_switch(&opt.skip_clip),
	 "skip reads that are soft/hard clipped")
	("paired_only,P", progopt::bool_switch(&opt.paired_only),
	 "output reads only if both matepair are chosen")
	("expand_read_id,E", progopt::bool_switch(&opt.expand_read_id),
	 "expand read identifiers by appending /1 or /2")
	("bam_file", progopt::value(&opt.bam_file)->required(), "BAM file as input");

	progopt::positional_options_description pos_desc;
	pos_desc.add("bam_file", 1);

	progopt::variables_map vm;
	try {
		progopt::store(progopt::command_line_parser(argc, argv)
				.options(desc)
				.positional(pos_desc).run(), vm);
		progopt::notify(vm);
	}
	catch(std::exception& e) {
		std::cerr << "ERROR: " << e.what() << std::endl;
		std::cerr << std::endl<< "Usage: " << argv[0] << " [options] <bam_file>" 
			<< std::endl << std::endl;
		std::cerr << desc << std::endl;
		exit(1);
	}

	if (vm.count("help")) {
		std::cout << std::endl<< "Usage: " << argv[0] << " [options] <bam_file>" 
			<< std::endl << std::endl;
		std::cout << desc << std::endl;
		exit(0);
	}

	if (vm.count("reference")) {
		opt.references = vm["reference"].as<std::vector<std::string>>();
	}

	if (!vm.count("output")) {
		std::string bam = vm["bam_file"].as<std::string>();
		opt.output_prefix = bam.substr(0, bam.find("."));
	}
}

int main(int argc, char *argv[])
{
	options opt;
	opt.skip_indel = false;
	opt.skip_nbase = false;
	opt.skip_clip = false;
	opt.paired_only = false;
	opt.expand_read_id = false;

	parse_options(argc, argv, opt);
	bam_filter(opt);
}
