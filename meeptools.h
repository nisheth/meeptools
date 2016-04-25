// include statements

#include <assert.h>
#include <limits.h>
#include <time.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "bseq.h"
#include "kseq.h"
#include "uthash.h"

// define statments

#define MEEPTOOLS_VERSION "0.2"

#define MAXREADLENGTH 500

#define STDSTRLEN 256

#define DEFAULTLCUT 35

#define MMODE_APPEND 0
#define MMODE_FILTER 1
#define MMODE_SORT 2
#define MMODE_STATS 3
#define MMODE_TRIM 4
#define MMODE_TRIM_SLOW 5
#define MMODE_TRIM_FAST 6

//  structures

typedef struct {
	unsigned long id;
	char *r1;
	char *r2;
	double mee1;
	double mee2;
	double rlen1;
	double rlen2;
	UT_hash_handle hh;
} readpair;

typedef struct {
	int chunk_size;
	int n_threads;
	double q2mee[75];
} opt_t;

typedef struct 
{
	unsigned long nreads;
	unsigned long nbases;
	double avgqsum;
	unsigned long qsumsum;
	double meesum;
	double meepsum;
	int minrl;
	int maxrl;
	float minaq;
	float maxaq;
	float minmeep;
	float maxmeep;
	unsigned long meep_bin_counts[7];
} read_stats;


typedef struct {
	const opt_t *opt;
	bseq_file_t *ks;
	bseq_file_t *ks2;
	read_stats in;
	read_stats in2;
	read_stats out;
	read_stats out2;
	read_stats singleout;
	int meeptools_mode;
	double mcut;
	int lcut;
	int paired; //1 is paired , 0 is single ended
	gzFile *outR1;
	gzFile *outR2;
	gzFile *outS;
	int gzipped;
	int truncate;
	time_t begin;
	unsigned long npairs;
	readpair *pairlist;
} shared_t;


typedef struct {
    int items;
	int *n_seqs; // number of reads or number of pairs
	int *n_seqs2;// number of reads or number of pairs
	bseq1_t **seqsarray; // read1
	bseq1_t **seqsarray2; // read2
	shared_t *es;
} step_t;

// function declarations

void bseq1_calculate_MEE_Qsum(bseq1_t *seq,double *q2mee);
void bseq1_update_MEE_Qsum(bseq1_t *seq,double *q2mee);
void bseq1_copy_MEE_Qsum(bseq1_t *seq);
int bseq1_MEE_trim(bseq1_t *seq,double *q2mee,int lcut,double mcut); // return 0 to filter out
int bseq1_MEE_trim_fast(bseq1_t *seq,double *q2mee,int lcut,double mcut); // return 0 to filter out
int bseq1_MEE_trim_slow(bseq1_t *seq,double *q2mee,int lcut,double mcut); // return 0 to filter out
char * char_repeat( int n, char c );
int copy_read(char *r,bseq1_t *s);
void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
void *kt_pipeline_cb(void *shared, int step, void *_data);
void read_stats_init(read_stats *rs);
int read_stats_update(read_stats *rs,bseq1_t *s,int when); // when=0 use o_ values else normal values... when=0 before, when=1 after
void read_stats_print(read_stats *rs);
int read_pair_comparator(readpair *rp, readpair *rp2);
void strnmcpy(char *s1, const char *s2, size_t m,size_t n);
static int subcommand_usage(int meeptools_mode);
static int usage();
int writeReadToGzFile(gzFile *gzf,bseq1_t *s,int gzipped);

