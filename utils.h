#ifndef COMMONBLOCK
#include "common.h"
#endif
#include<unistd.h>

int file_exists(const char* filename);

typedef struct
{
    unsigned long long int nreads;
    unsigned long long int nbases;
    double summee;
    double sumavgqual;
    unsigned int minRL;
    unsigned int maxRL;
    
    unsigned long long int nreads_meep1;
    unsigned long long int nreads_meep2;
    
    double avgRL;
    double overallMEEP;
    double avgRQ;
    double pmeep1; //percent reads with MEEP1
    double pmeep2; //percent reads with MEEP2
}readSetStats;

void init_q2mee_hash();
void seqIsInvalid(kseq_t *seq,char *fastqFilename);
int seqWriteToFileWithMateNumber(kseq_t *seq,gzFile fpout,double meep,double mee,double readQual,int write_mee,int write_readQual,int mate_number);
double seqCalculateMEE(kseq_t *seq);
double seqCalculateQScore(kseq_t *seq);
double seqCalculateQScoreExtra(kseq_t *seq,int *l70q20,int *l70q25,int *l70q30);
int readSetStatsInit(readSetStats *rss);
int readSetStatsAddRead(readSetStats *rss,int rl,double rq,double mee);
int readSetStatsUpdate(readSetStats *rss);
void readSetStatsPrintHeader();
void readSetStatsPrint(readSetStats *rss,readSetStats *rssBase,char *desc);
int str_split( char * str, char delim, char ***array, int *length );

