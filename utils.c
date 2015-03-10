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

int init_q2mee_hash(int offset)
{
    struct meeHash *s;
    int i;
    double mee;

    for (i=0;i<50;i++) {
        mee=pow(10.0,(-1.0*i/10));
        s=(struct meeHash*)malloc(sizeof(struct meeHash));
        s->Q=i+offset;
        s->MEE=mee;
        HASH_ADD_INT(q2mee,Q,s);
    }
    return 1;
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
    int j;
    double mee=0.0;
    struct meeHash *s;
    for (i=0;i<(int)strlen(seq->qual.s);i++) {
		j=(int)seq->qual.s[i];
        HASH_FIND_INT(q2mee,&j,s);
        mee+=s->MEE;
    }
    return mee;
}

double seqTrimCalculateMEEReadQuality(kseq_t *seq,int offset,int lcut,double mcut,int *newstart,int *newl,double *newReadQuality, double *untrimmedMEE, double *untrimmedReadQuality)
{
    int i,j,k;
    int readQual=0;
    int l=seq->qual.l;
    double mee[l];
    double summeeHash[l][l];
    double summeeCut[l];
    struct meeHash *s;
    double summee=0.0;
    int found=0;
    
    *newstart=-1;
    *newl=-1;
    
        
    for (i=0;i<l;i++) 
    {
        k=(int)seq->qual.s[i];
        HASH_FIND_INT(q2mee,&k,s);
        mee[i]=s->MEE;
		summee+=mee[i];
		readQual+=seq->qual.s[i];
		for (j=0;j<l;j++)
		{
			summeeHash[i][j]=-1.0;
		}
		summeeCut[i]=mcut*(l-i)/100;
    }
    
    *untrimmedMEE=summee;
    readQual-=(l*offset);
    *untrimmedReadQuality=readQual*1.0/l;
    
    if (summee < summeeCut[0])
    {
	*newstart=0;
	*newl=l;
	*newReadQuality=*untrimmedReadQuality;
	return summee;
    }

    summeeHash[0][0]=summee;
    summeeHash[1][0]=summee-mee[l-1];
    for (i=l-1;i>lcut-1;i--) 
    {
	for (j=0;j<l-i+1;j++) 
	{
	    if (j==0)
	    {
		summeeHash[l-(i-1)][j]=summeeHash[l-i][j]-mee[i-1];
	    }
	    else
	    {
	        summeeHash[l-i][j]=summeeHash[l-i][j-1]-mee[j-1]+mee[j+i-1];
	    }
			//printf("maxlen=%d,len=%d, l-i=%d, j=%d, summeHash=%.4f summeeCut=%.4f meep=%.4f mcut=%.4f, lcut=%d\n",maxlen,i,l-i,j,summeeHash[l-i][j],summeeCut[l-i],summeeHash[l-i][j]*100.0/i,mcut,lcut);
			//printf("len=%d, l-i=%d, j=%d, summeeHash=%.4f summeeCut=%.4f meep=%.4f mcut=%.4f\n",i,l-i,j,summeeHash[l-i][j],summeeCut[l-i],summeeHash[l-i][j]*100.0/i,mcut);
	    if ( summeeHash[l-i][j] < summeeCut[l-i] ) 
	    {
		found=1;
		*newstart=j;
		*newl=i;
	    }
	}
	if ( found == 1 ) 
        {
	    readQual=0;
	    for ( i = (*newstart) ; i <= (*newl) ; i++ )
	    {
	        readQual+=seq->qual.s[i];
	    }
	    readQual-=((*newl)*offset);
    	*newReadQuality=readQual*1.0/(*newl);
	    return summeeHash[l-*newl][*newstart];
	}
    }
    return -1.0;
}



double seqCalculateQScore(kseq_t *seq,int offset)
{
    unsigned int readQual=0,i,l;
    l=seq->qual.l;
    if (l==0) return 0.0;
    for (i=0;i<l;i++) {
        readQual+=seq->qual.s[i];
    }
    readQual-=(l*offset);
    return readQual*1.0/l;
}

double seqCalculateQScoreExtra(kseq_t *seq,int offset,int *l70q20,int *l70q25,int *l70q30)
{
    int readQual=0,i,l,q;
    int c20=0,c25=0,c30=0;
    l=seq->qual.l;
    if (l==0) return 0.0;
    for (i=0;i<l;i++) {
        q=seq->qual.s[i]-offset;
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
    fprintf(stdout,"\nNreads\tPercent_reads\tNbases\tPercent_bases\tminRL\tmaxRL\tavgRL\tavgRQ\toverallMEEP\t(Q equivalent)\tNreads_MEEP1\tNreads_MEEP2\tDescription\n\n");
}

void readSetStatsPrint(readSetStats *rss,readSetStats *rssBase,char *desc)
{
	fprintf(stdout,"%llu\t%.2f\t%llu\t%.2f\t%u\t%u\t%.2f\t%.2f\t%.4f\t%d\t%llu\t%llu\t%s\n", \
	rss->nreads, \
	(rss->nreads*1.0/rssBase->nreads)*100.0, \
	rss->nbases, \
	(rss->nbases*1.0/rssBase->nbases)*100.0, \
	rss->minRL, \
	rss->maxRL, \
	rss->avgRL, \
	rss->avgRQ, \
	rss->overallMEEP, \
	(int) (-10 * log10 (rss->overallMEEP/100.0)), \
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

int mykseqWriteToFile(mykseq *ks,gzFile *fpout,int write_mee,int write_readQual)
{
    char linebuffer[2 * MAXREADLENGTH];
    char commentMEEPtmp[STDSTRLEN];
    char commentMEEP[STDSTRLEN];
    
    
    if (write_mee && write_readQual) {
        sprintf(commentMEEPtmp,"MEEP=%.4f:MEE=%.4f:QUAL=%.2f",ks->meep1,ks->mee1,ks->rq1);
    }
    else if(write_mee)
    {
        sprintf(commentMEEPtmp,"MEEP=%.4f:MEE=%.4f",ks->meep1,ks->mee1);
    }
    else if(write_readQual)
    {
        sprintf(commentMEEPtmp,"MEEP=%.4f:QUAL=%.2f",ks->meep1,ks->rq1);
    }
    else
    {
        sprintf(commentMEEPtmp,"MEEP=%.4f",ks->meep1);
    }
    
    if ((int)strlen(ks->comment1)>0)
    {
        sprintf(commentMEEP,"%s:%s",ks->comment1,commentMEEPtmp);
    }
    else
    {
        sprintf(commentMEEP,"%s",commentMEEPtmp);
    }
     
    sprintf(linebuffer,"@%s %s\n%s\n+\n%s\n",ks->name1,commentMEEP,ks->bases1,ks->quality1);  
    gzwrite(fpout[0],linebuffer,strlen(linebuffer));
    
    if (fpout[1]!=NULL) 
    {
		if (write_mee && write_readQual) {
			sprintf(commentMEEPtmp,"MEEP=%.4f:MEE=%.4f:QUAL=%.2f",ks->meep2,ks->mee2,ks->rq2);
		}
		else if(write_mee)
		{
			sprintf(commentMEEPtmp,"MEEP=%.4f:MEE=%.4f",ks->meep2,ks->mee2);
		}
		else if(write_readQual)
		{
			sprintf(commentMEEPtmp,"MEEP=%.4f:QUAL=%.2f",ks->meep2,ks->rq2);
		}
		else
		{
			sprintf(commentMEEPtmp,"MEEP=%.4f",ks->meep2);
		}
	
		if ((int)strlen(ks->comment2)>0)
		{
			sprintf(commentMEEP,"%s:%s",ks->comment2,commentMEEPtmp);
		}
		else
		{
			sprintf(commentMEEP,"%s",commentMEEPtmp);
		}
	 
		sprintf(linebuffer,"@%s %s\n%s\n+\n%s\n",ks->name2,commentMEEP,ks->bases2,ks->quality2);  
		gzwrite(fpout[1],linebuffer,strlen(linebuffer));    
    }
    
    return 1;
}

int seqWriteSubseqToFileWithMateNumber(kseq_t *seq, int newstart, int newl, gzFile fpout, double meep, double mee, double readQual, int write_mee, int write_readQual, int mate_number)
{
    char linebuffer[2 * MAXREADLENGTH];
    char bases[MAXREADLENGTH];
    char quality[MAXREADLENGTH];
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
    
    
    memset(bases, '\0', sizeof(bases));
    memset(quality, '\0', sizeof(quality));
    
    strncpy(bases,&seq->seq.s[newstart],newl);
    strncpy(quality,&seq->qual.s[newstart],newl);
     
    sprintf(linebuffer,"@%s %s\n%s\n+\n%s\n",nametmp,commentMEEP,bases,quality);  
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



