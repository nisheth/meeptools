#include "meeptools_append.h"
static int meeptools_append_usage(int extra)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   meeptools append [options] \n\n Appends MEEP score to comment section of read id.\n\n");
    fprintf(stderr, "Options: -i FILE  FASTQ FILE(S) ... for paired-end read1 and read2 files should be comma separated\n");
    fprintf(stderr, "         -o FILE  FASTQ FILE(S) ... for paired-end read1 and read2 files should be comma separated\n");
    fprintf(stderr, "Optional:\n");
    fprintf(stderr, "         -f INT   Offset 33(default) or 64\n");
    fprintf(stderr, "         -q       Append read quality to read comment\n");
    fprintf(stderr, "         -m       Append MEE to read comment\n");
    fprintf(stderr, "         -h       help\n");
    fprintf(stderr, "\n");
    if (extra)
    {
        fprintf(stderr,"Examples:\n");
        fprintf(stderr,"meeptools append -i a.fastq.gz -o a_meep.fastq.gz\n");
        fprintf(stderr,"meeptools append -i read1.fastq,read2.fastq -o read1_meep.fastq.gz,read2_meep.fastq.gz -f 64\n");
    }
    return 1;
}

int meeptools_append(int argc, char *argv[])
{
    int c,i,l,offset=33;
    int iflag=0;
    int oflag=0;
    int qflag=0;
    int mflag=0;
    char *inputFastqstr;
    char **inputFastqs;
    int nInputFastqs;
    char *outputFastqstr;
    char **outputFastqs;
    int nOutputFastqs;
    gzFile fp,fpout;
    kseq_t *seq;
    double mee=0.0;
    double readQual=0.0;
    double meep=0.0;
    int nreads;
    
    char msg[STDSTRLEN];
    
    while ((c = getopt(argc, argv, "i:f:o:mqh")) != -1)
    {
        switch (c)
        {
        case 'i':
            iflag = 1;
            inputFastqstr = strdup(optarg);
            break;
        case 'f':
        	offset = atoi(optarg);
        	break;
        case 'o':
            oflag = 1;
            outputFastqstr = strdup(optarg);
            break;
        case 'q':
            qflag = 1;
            break;
        case 'm':
            mflag = 1;
            break;
        case 'h':
            return meeptools_append_usage(1);
        default:
            return meeptools_append_usage(0);
        }
    }
    if (offset !=33 && offset !=64 ) {
    	ErrorMsgExit("Offset value can only be 33(default) or 64!");
    }
    if (!init_q2mee_hash(offset))
    {
        ErrorMsgExit("Initialization of q2mee failed!");
    }
	DebugMsg("Initializing q2mee hash complete.");
    if (iflag==0)
    {
        ErrorMsg("missing -i option");
        return meeptools_append_usage(0);
    }
    if (oflag==0)
    {
        ErrorMsg("missing -o option");
        return meeptools_append_usage(0);
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
    
    sprintf(msg,"Number of input Fastq(s) : %d",nInputFastqs);
    DebugMsg(msg);
    
    for (i=0;i<nInputFastqs;i++)
    {
        if ((fp = gzopen(inputFastqs[i], "rb")) == NULL )
        {
            sprintf(msg, "%s file cannot be opened!",inputFastqs[i]);
            ErrorMsgExit(msg);    
        } 
        seq = kseq_init(fp);
        nreads=0;
        sprintf(msg,"FASTQ file %s opened.",inputFastqs[i]);
        DebugMsg(msg);
        
        if ((fpout = gzopen(outputFastqs[i], "wb")) == NULL )
        {
            sprintf(msg, "%s file cannot be opened!",outputFastqs[i]);
            ErrorMsgExit(msg);    
        } 
        
        sprintf(msg,"FASTQ file %s opened.",outputFastqs[i]);
        DebugMsg(msg);
        
        while (0==0)
        {
            l = kseq_read(seq);
            if (l == -1) break;
            if (l == -2)
            {
                ErrorMsg("Invalid sequence detected!");
                seqIsInvalid(seq,inputFastqs[i]);
            }

            mee = seqCalculateMEE(seq);
            meep = mee*100.0/l;
            if (qflag) readQual = seqCalculateQScore(seq,offset);
            
            if (!seqWriteToFileWithMateNumber(seq,fpout,meep,mee,readQual,mflag,qflag,0)) {
                sprintf(msg,"Writing seqid %s to file %s failed!",seq->name.s,outputFastqs[i]);
                ErrorMsgExit(msg);
            }
            
            if(nreads%PRINTINTERVAL == 0)
            {
                sprintf(msg,"%d reads processed.",nreads);
                PrintMsg(msg);
            }
            nreads++;
        
        }
        
        kseq_destroy(seq);
        
        gzclose(fp);
        
        sprintf(msg,"FASTQ file %s closed.",inputFastqs[i]);
        DebugMsg(msg);
        
        gzclose(fpout);
        
        sprintf(msg,"FASTQ file %s closed.",outputFastqs[i]);
        DebugMsg(msg);
    
    }
    
    
    return 1;
}
