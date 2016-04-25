#ifndef BFC_BSEQ_H
#define BFC_BSEQ_H

// includes

#include <stdint.h>

// defines

#define NO_OUTPUT 0 // do not output this read or pair
#define PAIRED_OUTPUT 1
#define SINGLE_OUTPUT 2
// define output modes... valid only for paired end input
// pairoutput means read1/2 input read goes to read1/2 output
// singleoutput means reads1/2 goes to single ended output file if provided... ie. one of the pair does not satisfy the filter or trim criteria

// structs

struct bseq_file_s;

typedef struct bseq_file_s bseq_file_t;

typedef struct {
	int l_seq;
	char *name, *comment, *seq, *qual;
	double o_mee,o_meep; // original before trimming values are all "o_"
	double mee,meep;
	int qsum;
	int o_qsum;
	double avgq;
	double o_avgq;
	int trim_start;
	int trim_len;
	int output_mode;
} bseq1_t;

//  functions

bseq_file_t *bseq_open(const char *fn);
void bseq_close(bseq_file_t *fp);
bseq1_t *bseq_read(bseq_file_t *fp, int chunk_size, int keep_comment, int *n_);

#endif
