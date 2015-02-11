#include "meeptools_filter.h"
static int meeptools_filter_usage(int extra)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   meeptools filter [options] \n\n Filters reads by MEEP score.\n\n");
    fprintf(stderr, "Options: -f FILE         FASTQ FILE(S) ... for paired-end read1 and read2 files should be comma separated\n");
    fprintf(stderr, "         -o FILE         FASTQ FILE(S) ... for paired-end read1 and read2 files should be comma separated\n");
    fprintf(stderr, "         -m FLOAT        MEEP score cut off (between 0 and 20)\n");
    fprintf(stderr, "         -l INTEGER      read length cut off (optional)\n");
    fprintf(stderr, "         -s FILE         FASTQ FILE for single read (read1 or read2) meeting MEEP score cut off\n");
    fprintf(stderr, "         -h       help\n");
    fprintf(stderr, "\n");
    if (extra)
    {
        fprintf(stderr,"Examples:\n");
        fprintf(stderr,"meeptools filter -f a.fastq.gz -o a_mfiltered.fastq.gz -m 1.0\n");
        fprintf(stderr,"meeptools filter -f read1.fastq,read2.fastq -o read1_meep_filtered.fastq.gz,read2_meep_filtered.fastq.gz -m 1.0\n");
        fprintf(stderr,"meeptools filter -m 1.0 -l 50 -f read1.fastq,read2.fastq -o read1_meep_filtered.fastq.gz,read2_meep_filtered.fastq.gz -s single_end_meep_filtered.fastq.gz\n");
    }
    return 1;
}
int meeptools_filter(int argc, char *argv[])
{
    int c,i;
    int l[2];
    int fflag=0;
    int oflag=0;
    int mflag=0;
    int sflag=0;
    int lcut=0;
    double mcut=0.0;
    char *inputFastqstr;
    char **inputFastqs;
    int nInputFastqs;
    char *outputFastqstr;
    char **outputFastqs;
    char *singleEndFastqOutFile;
    int nOutputFastqs;
    gzFile fp[2],fpout[2];
    gzFile sfp;
    kseq_t *seq1,*seq2;
    double mee[2];
    double meep[2];
    int nreads;
    
    char msg[STDSTRLEN];
    
    for (i=0;i<2;i++) {
        mee[i]=0.0;
        meep[i]=0.0;
    }
    
    while ((c = getopt(argc, argv, "f:o:m:l:s:h")) != -1)
    {
        switch (c)
        {
        case 'f':
            fflag = 1;
            inputFastqstr = strdup(optarg);
            break;
        case 'o':
            oflag = 1;
            outputFastqstr = strdup(optarg);
            break;
        case 's':
            sflag = 1;
            singleEndFastqOutFile = strdup(optarg);
            break;
        case 'm':
            mflag=1;
            mcut = strtod(optarg,&msg);
            break;
        case 'l':
            lcut = atoi(optarg);
            break;
        case 'h':
            return meeptools_filter_usage(1);
        default:
            return meeptools_filter_usage(0);
        }
    }
    
    sprintf(msg,"MEEP cut off = %.2f",mcut);
    DebugMsg(msg);
    sprintf(msg,"RL cut off = %d",lcut);
    DebugMsg(msg);
    
    if (fflag==0)
    {
        ErrorMsg("missing -f option");
        return meeptools_filter_usage(0);
    }
    if (oflag==0)
    {
        ErrorMsg("missing -o option");
        return meeptools_filter_usage(0);
    }
    if (lcut<0)
    {
        ErrorMsg("Invalid -l option");
        return meeptools_filter_usage(0);
    }
    if (mflag==0)
    {
        ErrorMsg("missing -m option");
        return meeptools_filter_usage(0);
    }
    else if (mcut<0 || mcut>20)
    {
        ErrorMsg("MEEP cut off value out of bounds!");
        return meeptools_filter_usage(0);
    }
    if (!str_split( inputFastqstr, ',', &inputFastqs, &nInputFastqs ))
    {
        ErrorMsgExit("Spliting of delimited input failed!");
    }
    if (!str_split( outputFastqstr, ',', &outputFastqs, &nOutputFastqs ))
    {
        ErrorMsgExit("Spliting of delimited output failed!");
    }
    if (nInputFastqs!=nOutputFastqs)
    {
        ErrorMsg("Number of Input and Output Fastq(s) do not match!");
        sprintf(msg,"Number of Input Fastq(s) = %d, Number of Output Fastq(s) = %d.",nInputFastqs,nOutputFastqs);
        ErrorMsgExit(msg);
    }
    if (nInputFastqs>2)
    {
        ErrorMsg("Number of Input Fastq(s) cannot exceed 2!");
        sprintf(msg,"Number of Input Fastq(s) = %d, Number of Output Fastq(s) = %d.",nInputFastqs,nOutputFastqs);
        ErrorMsgExit(msg);        
    }
    if (sflag && nInputFastqs==1) {
        ErrorMsg("-s options invalid with just one input fastq.");
        return meeptools_filter_usage(1);
    }
    
    if ((fp[0] = gzopen(inputFastqs[0], "rb")) == NULL )
    {
        sprintf(msg, "%s file cannot be opened!",inputFastqs[0]);
        ErrorMsgExit(msg);    
    }
    seq1 = kseq_init(fp[0]);
    nreads=0;
    
    sprintf(msg,"FASTQ file %s opened.",inputFastqs[0]);
    DebugMsg(msg);
    
    if ((fpout[0] = gzopen(outputFastqs[0], "wb")) == NULL )
    {
        sprintf(msg, "%s file cannot be opened!",outputFastqs[0]);
        ErrorMsgExit(msg);    
    } 
    
    sprintf(msg,"FASTQ file %s opened.",outputFastqs[0]);
    DebugMsg(msg);
    
    if (nInputFastqs==2)
    {
        if ((fp[1] = gzopen(inputFastqs[1], "rb")) == NULL )
        {
            sprintf(msg, "%s file cannot be opened!",inputFastqs[1]);
            ErrorMsgExit(msg);    
        }
        seq2 = kseq_init(fp[1]);        
    
        sprintf(msg,"FASTQ file %s opened.",inputFastqs[1]);
        DebugMsg(msg);
        
        if ((fpout[1] = gzopen(outputFastqs[1], "wb")) == NULL )
        {
            sprintf(msg, "%s file cannot be opened!",outputFastqs[1]);
            ErrorMsgExit(msg);    
        } 
        
        sprintf(msg,"FASTQ file %s opened.",outputFastqs[1]);
        DebugMsg(msg);
        
        if(sflag)
        {
            if ((sfp = gzopen(singleEndFastqOutFile, "wb")) == NULL )
            {
                sprintf(msg, "%s file cannot be opened!",singleEndFastqOutFile);
                ErrorMsgExit(msg);    
            } 
            
            sprintf(msg,"FASTQ file %s opened.",singleEndFastqOutFile);
            DebugMsg(msg);            
        }
    }
     
    while (0==0)
    {
        l[0] = kseq_read(seq1);
        if (l[0] == -1) break;
        if (l[0] == -2)
        {
            ErrorMsg("Invalid sequence detected!");
            seq_is_invalid(seq1,inputFastqs[0]);
        }
        mee[0] = calculate_mee(seq1);
        meep[0] = mee[0]*100.0/l[0];
        if (nInputFastqs==2)
        {
            l[1] = kseq_read(seq2);
            if (l[1] == -1) break;
            if (l[1] == -2)
            {
                ErrorMsg("Invalid sequence detected!");
                seq_is_invalid(seq2,inputFastqs[1]);
            }            
            mee[1] = calculate_mee(seq2);
            meep[1] = mee[1]*100.0/l[1];
        }
        if ((meep[0] < mcut) && (meep[1] < mcut))
        {
            if (!seq_write_to_file(seq1,fpout[0],meep[0],0,0,0,0)) {
                sprintf(msg,"Writing seqid %s to file %s failed!",seq1->name.s,outputFastqs[0]);
                ErrorMsgExit(msg);
            }
            if (!seq_write_to_file(seq2,fpout[1],meep[1],0,0,0,0)) {
                sprintf(msg,"Writing seqid %s to file %s failed!",seq2->name.s,outputFastqs[1]);
                ErrorMsgExit(msg);
            }
            continue;
        }
        if ((meep[0] < mcut) && sflag)
        {
            if (!seq_write_to_file(seq1,sfp,meep[0],0,0,0,0)) {
                sprintf(msg,"Writing seqid %s to file %s failed!",seq1->name.s,singleEndFastqOutFile);
                ErrorMsgExit(msg);
            }
            continue;
        }
        if ((meep[1] < mcut) && sflag)
        {
            if (!seq_write_to_file(seq2,sfp,meep[0],0,0,0,0)) {
                sprintf(msg,"Writing seqid %s to file %s failed!",seq2->name.s,singleEndFastqOutFile);
                ErrorMsgExit(msg);
            }
            continue;
        }
        if ((meep[0] < mcut) && !sflag)
        {
            if (!seq_write_to_file(seq1,fpout[0],meep[0],0,0,0,0)) {
                sprintf(msg,"Writing seqid %s to file %s failed!",seq1->name.s,outputFastqs[0]);
                ErrorMsgExit(msg);
            }
            continue;
        }
    }
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    gzclose(fp[0]);
    gzclose(fpout[0]);
    if (nInputFastqs==2) {
        gzclose(fp[1]);
        gzclose(fpout[1]);
        if (sflag) gzclose(sfp);
    }
    DebugMsg("Done closing all files");
    return 1;
}
