#ifndef _BAM_SUBSAMPLE_H_
#define _BAM_SUBSAMPLE_H_

#include <vector>
#include <string>

// convert char type ACTG to 2-bit representation:
// A -> 0, C -> 1, T -> 2, G -> 3
#define base_to_bits(b) (((b) >> 1) & 3)

// seperator used to expand read identifier
#define READ_EXPAND_SEP "/"

// complementary bases in PREMIER's order
const char pmr_complements[4] = {'T', 'G', 'A', 'C'};
// nucleotides in samtools's order
const char sam_nucleotides[5] = {'A', 'C', 'G', 'T', 'N'};

const char aux_md[2] = {'M', 'D'};

typedef struct _options {
	bool skip_indel;
	bool skip_clip;
	bool skip_nbase;
	bool paired_only;
	bool expand_read_id;

	std::vector<std::string> references;
	int cutoff_start;
	int cufoff_end;

	std::string output_prefix;
	std::string bam_file;
} options;

#endif
