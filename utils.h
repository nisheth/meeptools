#ifndef COMMONBLOCK
#include "common.h"
#endif
#include<unistd.h>
//#include <stdlib.h>

#ifndef NBINS 31
#define NBINS 31
#endif

int file_exists(const char* filename);

typedef struct
{
    unsigned long long int nreads;
    unsigned long long int nbases;
    double summee;
    double sumavgqual;
    unsigned int minRL;
    unsigned int maxRL;
    unsigned int meepbins[31]; // 31 bins with binwidth of 0.1 meep    
    unsigned long long int nreads_meep1;
    unsigned long long int nreads_meep2;
    
    double avgRL;
    double overallMEEP;
    double avgRQ;
    double pmeep1; //percent reads with MEEP1
    double pmeep2; //percent reads with MEEP2
}readSetStats;

typedef struct 
{
	int count; //key
	char *name1;
	char *name2;
	char *comment1;
	char *comment2;
	char *bases1;
	char *bases2;
	char *quality1;
	char *quality2;
	double mee1;
	double mee2;
	int l1;
	int l2;
	double rq1;
	double rq2;
	double meep1;
	double meep2;
	double combinedMEEP;
	UT_hash_handle hh;
}mykseq;

//int floatcomp(const void* elem1, const void* elem2);
int init_q2mee_hash(int offset);
void seqIsInvalid(kseq_t *seq,char *fastqFilename);
int seqWriteToFileWithMateNumber(kseq_t *seq,gzFile fpout,double meep,double mee,double readQual,int write_mee,int write_readQual,int mate_number);
int seqWriteSubseqToFileWithMateNumber(kseq_t *seq, int newstart, int newl, gzFile fpout, double meep, double mee, double readQual, int write_mee, int write_readQual, int mate_number);
double seqCalculateMEE(kseq_t *seq);
double seqCalculateQScore(kseq_t *seq,int offset);
double seqCalculateQScoreExtra(kseq_t *seq,int offset,int *l70q20,int *l70q25,int *l70q30);
double seqTrimCalculateMEEReadQuality(kseq_t *seq,int offset,int lcut,double mcut,int *newstart,int *newl,double *newReadQuality, double *untrimmedlMEE, double *untrimmedReadQuality);
int readSetStatsInit(readSetStats *rss);
int readSetStatsAddRead(readSetStats *rss,int rl,double rq,double mee);
int readSetStatsUpdate(readSetStats *rss);
void readSetStatsPrintHeader();
void readSetStatsPrintHist(readSetStats *rss,char *fastqFilename);
void readSetStatsPrintJson(readSetStats *rss,char *fastqFilename);
void readSetStatsPrint(readSetStats *rss,readSetStats *rssBase,char *desc);
int str_split( char * str, char delim, char ***array, int *length );
int mykseqWriteToFile(mykseq *ks,gzFile *fpout,int write_mee,int write_readQual);

