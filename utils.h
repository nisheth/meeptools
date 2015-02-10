#ifndef COMMONBLOCK
#include "common.h"
#endif
#include<unistd.h>

int file_exists(const char* filename);

void init_q2mee_hash();
void seq_is_invalid(kseq_t *seq,char *fastqFilename);
double calculate_mee(kseq_t *seq);
double calculate_qscore(kseq_t *seq);
double calculate_qscore_extra(kseq_t *seq,int *l70q20,int *l70q25,int *l70q30);
