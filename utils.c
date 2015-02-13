#include "utils.h"
extern struct meeHash *q2mee=NULL;

int file_exists(const char* filename)
{
    int result;
    result = access (filename, F_OK); // F_OK tests existence also (R_OK,W_OK,X_OK).
                                      //            for readable, writeable, executable
    if ( result == 0 )
    {
       return 1;
    }
    else
    {
       return 0;
    }
}

void init_q2mee_hash()
{
    struct meeHash *s;
    char tmpq;
    char q[2];
    int i;
    double mee;

    for (i=0;i<50;i++) {
        tmpq=(char)(i+33);
        sprintf(q,"%c\0",tmpq);
        mee=pow(10.0,(-1.0*i/10));
        s=(struct meeHash*)malloc(sizeof(struct meeHash));
        strcpy(s->Q,q);
        s->MEE=mee;
        HASH_ADD_STR(q2mee,Q,s);
    }
}

void seqIsInvalid(kseq_t *seq,char *fastqFilename)
{
    char msg[STDSTRLEN];
    sprintf(msg,"seqid %s, sequence string length=%d, quality string length=%d.",seq->name.s,(int)strlen(seq->seq.s),(int)strlen(seq->qual.s));
    ErrorMsgExit(msg);
}

double seqCalculateMEE(kseq_t *seq)
{
    int i;
    char q[2];
    double mee=0.0;
    struct meeHash *s;
    for (i=0;i<(int)strlen(seq->qual.s);i++) {
        sprintf(q,"%c\0",seq->qual.s[i]);
        HASH_FIND_STR(q2mee,q,s);
        mee+=s->MEE;
    }
    return mee;
}

double seqCalculateQScore(kseq_t *seq)
{
    unsigned int readQual=0,i,l;
    l=seq->qual.l;
    if (l==0) return 0.0;
    for (i=0;i<l;i++) {
//        readQual+=(seq->qual.s[i]-'!');
        readQual+=seq->qual.s[i];
    }
    readQual-=(l*'!');
    return readQual*1.0/l;
}

double seqCalculateQScoreExtra(kseq_t *seq,int *l70q20,int *l70q25,int *l70q30)
{
    int readQual=0,i,l,q;
    int c20=0,c25=0,c30=0;
    l=seq->qual.l;
    if (l==0) return 0.0;
    for (i=0;i<l;i++) {
        q=seq->qual.s[i]-'!';
        readQual+=q;
        if (q>=20) c20++;
        if (q>=25) c25++;
        if (q>=30) c30++; 
    }
    if (c20/l>=0.7) {*l70q20=1;} else {*l70q20=0;}
    if (c25/l>=0.7) {*l70q25=1;} else {*l70q25=0;}
    if (c30/l>=0.7) {*l70q30=1;} else {*l70q30=0;}
    return readQual/l;
}

int readSetStatsInit(readSetStats *rss)
{
    rss->nreads=0;
    rss->nbases=0;
    rss->summee=0.0;
    rss->sumavgqual=0.0;
    rss->minRL=(unsigned int)MAXREADLENGTH;
    rss->maxRL=0;
    
    rss->avgRL=0.0;
    rss->overallMEEP=0.0;
    rss->avgRQ=0.0;
    rss->nreads_meep1=0;
    rss->nreads_meep2=0;
    rss->pmeep1=0.0;
    rss->pmeep2=0.0;
    
    return 1;
}

int readSetStatsAddRead(readSetStats *rss,int rl,double rq,double mee)
{
    rss->nreads++;
    rss->nbases+=rl;
    rss->summee+=mee;
    rss->sumavgqual+=rq;
    
    if (rl>rss->maxRL) rss->maxRL=rl;
    if (rl<rss->minRL) rss->minRL=rl;
    if (mee*100.0/rl<1) rss->nreads_meep1++;
    if (mee*100.0/rl<2) rss->nreads_meep2++;
    
    return 1;
}

int readSetStatsUpdate(readSetStats *rss)
{
    if (rss->nreads>0)
    {
        rss->avgRL=rss->nbases*1.0/rss->nreads;
        rss->avgRQ=rss->sumavgqual/rss->nreads;
        rss->pmeep1=(rss->nreads_meep1*1.0/rss->nreads)*100.0;
        rss->pmeep2=(rss->nreads_meep2*1.0/rss->nreads)*100.0;     
    }
    if (rss->nbases>0) rss->overallMEEP=rss->summee*100.0/rss->nbases;
    
    return 1;
}

void readSetStatsPrintHeader()
{
    fprintf(stdout,"\nNreads\tPercent_reads\tNbases\tPercent_bases\tminRL\tmaxRL\tavgRL\tavgRQ\toverallMEEP\tNreads_MEEP1\tNreads_MEEP2\tDescription\n\n");
}

void readSetStatsPrint(readSetStats *rss,readSetStats *rssBase,char *desc)
{
	fprintf(stdout,"%llu\t%.2f\t%llu\t%.2f\t%u\t%u\t%.2f\t%.2f\t%.4f\t%llu\t%llu\t%s\n", \
	rss->nreads, \
	(rss->nreads*1.0/rssBase->nreads)*100.0, \
	rss->nbases, \
	(rss->nbases*1.0/rssBase->nbases)*100.0, \
	rss->minRL, \
	rss->maxRL, \
	rss->avgRL, \
	rss->avgRQ, \
	rss->overallMEEP, \
	rss->nreads_meep1, \
	rss->nreads_meep2, \
	desc);
}


int seqWriteToFileWithMateNumber(kseq_t *seq,gzFile fpout,double meep,double mee,double readQual,int write_mee,int write_readQual,int mate_number)
{
    char linebuffer[2 * MAXREADLENGTH];
    char commentMEEPtmp[STDSTRLEN];
    char commentMEEP[STDSTRLEN];
    char nametmp[STDSTRLEN];
    char word[STDSTRLEN];
    
    if (mate_number < 0 || mate_number >2) ErrorMsgExit("Invalid mate number supplied!");
    
    if (mate_number==1 || mate_number==2)
    {
	sprintf(word,"/%d",mate_number);
	if(strstr(seq->name.s,word)==NULL) sprintf(nametmp,"%s/%d",seq->name.s,mate_number);
    }
    else
    {
	sprintf(nametmp,"%s",seq->name.s);
    }
    
    if (write_mee && write_readQual) {
        sprintf(commentMEEPtmp,"MEEP=%.4f:MEE=%.4f:QUAL=%.2f",meep,mee,readQual);
    }
    else if(write_mee)
    {
        sprintf(commentMEEPtmp,"MEEP=%.4f:MEE=%.4f",meep,mee);
    }
    else if(write_readQual)
    {
        sprintf(commentMEEPtmp,"MEEP=%.4f:QUAL=%.2f",meep,readQual);
    }
    else
    {
        sprintf(commentMEEPtmp,"MEEP=%.4f",meep);
    }
    
    if (seq->comment.l)
    {
        sprintf(commentMEEP,"%s:%s",seq->comment.s,commentMEEPtmp);
    }
    else
    {
        sprintf(commentMEEP,"%s",commentMEEPtmp);
    }
     
    sprintf(linebuffer,"@%s %s\n%s\n+\n%s\n",nametmp,commentMEEP,seq->seq.s,seq->qual.s);  
    gzwrite(fpout,linebuffer,strlen(linebuffer));
    
    return 1;
}


int str_split( char * str, char delim, char ***array, int *length ) {
  char *p;
  char **res;
  int count=0;
  int k=0;

  p = str;
  // Count occurance of delim in string
  while( (p=strchr(p,delim)) != NULL ) {
    *p = 0; // Null terminate the deliminator.
    p++; // Skip past our new null
    count++;
  }

  // allocate dynamic array
  res = calloc( 1, (count+1) * sizeof(char *));
  if( !res ) return -1;

  p = str;
  for( k=0; k<count; k++ ){
    if( *p ) res[k] = p;  // Copy start of string
    p = strchr(p, 0 );    // Look for next null
    p++; // Start of next string
  }
  res[count]=p;
  *array = res;
  *length = count+1;

  return 1;
}
/*
example:
char str[] = "JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC,";

int main() {
  char **res;
  int k=0;
  int count =0;
  int rc;

  rc = str_split( str, ',', &res, &count );
  if( rc ) {
    printf("Error: %s errno: %d \n", strerror(errno), errno);
  }

  printf("count: %d\n", count );
  for( k=0; k<count; k++ ) {
    printf("str: %s\n", res[k]);
  }

  free(res );
  return 0;
}
 */


