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

void seq_is_invalid(kseq_t *seq,char *fastqFilename)
{
    fprintf(stderr,"[%s]: seqid %s, sequence string length=%d, quality string length=%d\n",__func__,seq->name.s,(int)strlen(seq->seq.s),(int)strlen(seq->qual.s));
    exit(1);
}

double calculate_mee(kseq_t *seq)
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

double calculate_qscore(kseq_t *seq)
{
    int readQual=0,i,l;
    l=seq->qual.l;
    if (l==0) return 0.0;
    for (i=0;i<l;i++) {
        readQual+=(seq->qual.s[i]-'!');
    }
    return readQual/l;
}

double calculate_qscore_extra(kseq_t *seq,int *l70q20,int *l70q25,int *l70q30)
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
/*
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
  res = calloc( 1, count * sizeof(char *));
  if( !res ) return -1;

  p = str;
  for( k=0; k<count; k++ ){
    if( *p ) res[k] = p;  // Copy start of string
    p = strchr(p, 0 );    // Look for next null
    p++; // Start of next string
  }

  *array = res;
  *length = count;

  return 0;
}

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

