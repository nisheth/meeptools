#include "meeptools_trim.h"
static int meeptools_trim_usage(int extra)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   meeptools trim [options] \n\n Trims reads from 5' and 3' ends to meet MEEP cutoff.\n\n");
    fprintf(stderr, "Options: -f FILE         FASTQ FILE(S) ... for paired-end read1 and read2 files should be comma separated\n");
    fprintf(stderr, "         -o FILE         FASTQ FILE(S) ... for paired-end read1 and read2 files should be comma separated\n");
    fprintf(stderr, "         -c FLOAT        MEEP score cut off (between 0 and 20)\n");
    fprintf(stderr, "Optional:\n");
    fprintf(stderr, "         -l INTEGER      read length cut off (default 35)\n");
    fprintf(stderr, "         -s FILE         FASTQ FILE for single read (read1 or read2) meeting MEEP score cut off\n");
    fprintf(stderr, "         -q              Append read quality to read comment\n");
    fprintf(stderr, "         -m              Append MEE to read comment\n");
    fprintf(stderr, "         -h              help\n");
    fprintf(stderr, "\n");
    if (extra)
    {
        fprintf(stderr,"Examples:\n");
        fprintf(stderr,"meeptools trim -f a.fastq.gz -o a_mtrimed.fastq.gz -m 1.0\n");
        fprintf(stderr,"meeptools trim -f read1.fastq,read2.fastq -o read1_meep_trimed.fastq.gz,read2_meep_trimed.fastq.gz -m 1.0\n");
        fprintf(stderr,"meeptools trim -m 1.0 -l 50 -f read1.fastq,read2.fastq -o read1_meep_trimed.fastq.gz,read2_meep_trimed.fastq.gz -s single_end_meep_trimed.fastq.gz\n");
    }
    return 1;
}
int meeptools_trim(int argc, char *argv[])
{
    int z,i;
    int l[2];
    int fflag=0;
    int oflag=0;
    int cflag=0;
    int mflag=0;
    int qflag=0;
    int sflag=0;
    int lcut=DEFAULTLCUT;
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
    double readQual[2];
    readSetStats rssIn[2];
    readSetStats rssOut[3];
    int nreads;
    int newstart[2],newl[2];
    int newstarttmp,newltmp;
    double newReadQuality,untrimmedMEE,untrimmedReadQuality;
    
    char msg[STDSTRLEN];
    
    for (i=0;i<2;i++) 
    {
        mee[i]=-1.0;
        meep[i]=-1.0;
        readQual[i]=0.0;
        newstart[i]=-1;
        newl[i]=-1;
        if (!readSetStatsInit(&rssIn[i])) ErrorMsgExit("failed initialization of read set!");
        if (!readSetStatsInit(&rssOut[i])) ErrorMsgExit("failed initialization of read set!");
    }
    
    while ((z = getopt(argc, argv, "f:o:c:l:s:mqh")) != -1)
    {
        switch (z)
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
        case 'c':
            cflag=1;
            mcut = strtod(optarg,&msg);
            break;
        case 'm':
            mflag=1;
            break;
        case 'q':
            qflag=1;
            break;
	case 'l':
            lcut = atoi(optarg);
            break;
        case 'h':
            return meeptools_trim_usage(1);
        default:
            return meeptools_trim_usage(0);
        }
    }

    sprintf(msg,"MEEP cut off = %.2f",mcut);
    DebugMsg(msg);
    sprintf(msg,"RL cut off = %d",lcut);
    DebugMsg(msg);
    
    if (fflag==0)
    {
        ErrorMsg("missing -f option");
        return meeptools_trim_usage(0);
    }
    if (oflag==0)
    {
        ErrorMsg("missing -o option");
        return meeptools_trim_usage(0);
    }
    if (lcut<0)
    {
        ErrorMsg("Invalid -l option");
        return meeptools_trim_usage(0);
    }
    if (cflag==0)
    {
        ErrorMsg("missing -c option");
        return meeptools_trim_usage(0);
    }
    else if (mcut<0 || mcut>20)
    {
        ErrorMsg("MEEP cut off value out of bounds!");
        return meeptools_trim_usage(0);
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
        return meeptools_trim_usage(1);
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

			if (!readSetStatsInit(&rssOut[2])) ErrorMsgExit("failed initialization of read set!");
            
            sprintf(msg,"FASTQ file %s opened.",singleEndFastqOutFile);
            DebugMsg(msg);            
        }
    }
    
    while (0==0)
    {
        for (i=0;i<2;i++) 
    	{
			mee[i]=-1.0;
			readQual[i]=0.0;
			newstart[i]=-1;
			newl[i]=-1;
        }
        newstarttmp=-1;
        newltmp=-1;
        
        l[0] = kseq_read(seq1);
        if (l[0] == -1) break;
        if (l[0] == -2)
        {
            ErrorMsg("Invalid sequence detected!");
            seqIsInvalid(seq1,inputFastqs[0]);
        }
        
        mee[0]=seqTrimCalculateMEEReadQuality(seq1,lcut,mcut,&newstarttmp,&newltmp,&newReadQuality,&untrimmedMEE,&untrimmedReadQuality);
        readQual[0]=newReadQuality;
        newstart[0]=newstarttmp;
        newl[0]=newltmp;
        meep[0]=mee[0]*100.0/newl[0];
        
        if (!readSetStatsAddRead(&rssIn[0],l[0],untrimmedReadQuality,untrimmedMEE)) ErrorMsgExit("failed adding read to read set!");
        if (nInputFastqs==2)
        {
            l[1] = kseq_read(seq2);
            if (l[1] == -1) break;
            if (l[1] == -2)
            {
                ErrorMsg("Invalid sequence detected!");
                seqIsInvalid(seq2,inputFastqs[1]);
            }            
            mee[1] = seqTrimCalculateMEEReadQuality(seq2,lcut,mcut,&newstarttmp,&newltmp,&newReadQuality,&untrimmedMEE,&untrimmedReadQuality);
	   		readQual[1] = newReadQuality;
	   		newstart[1] = newstarttmp;
	   		newl[1] = newltmp;
	   		meep[1]=mee[1]*100.0/newl[1];
	   		
	    	if (!readSetStatsAddRead(&rssIn[1],l[1],untrimmedReadQuality,untrimmedMEE)) ErrorMsgExit("failed adding read to read set!");
        }
        
        if(nreads%PRINTINTERVAL == 0)
		{
			sprintf(msg,"%d reads processed.",nreads);
			PrintMsg(msg);
		}
		nreads++;
        
		if (mee[0]!=-1 && mee[1]!=-1) 
		{
        	if (!readSetStatsAddRead(&rssOut[0],newl[0],readQual[0],mee[0])) ErrorMsgExit("failed adding read to read set!");
        	if (!readSetStatsAddRead(&rssOut[1],newl[1],readQual[1],mee[1])) ErrorMsgExit("failed adding read to read set!");
			continue;
		}
		if (mee[0]!=-1 && !sflag)
		{
		    if (!seqWriteSubseqToFileWithMateNumber(seq1,newstart[0],newl[0],fpout[0],meep[0],mee[0],readQual[0],mflag,qflag,0)) {
                sprintf(msg,"Writing seqid %s to file %s failed!",seq1->name.s,outputFastqs[0]);
                ErrorMsgExit(msg);
            }
        	if (!readSetStatsAddRead(&rssOut[0],newl[0],readQual[0],mee[0])) ErrorMsgExit("failed adding read to read set!");
			continue;
		}
		if (mee[0]!=-1 && sflag)
		{
        	if (!readSetStatsAddRead(&rssOut[2],newl[0],readQual[0],mee[0])) ErrorMsgExit("failed adding read to read set!");
			continue;
		}
		if (mee[1]!=-1 && sflag)
		{
        	if (!readSetStatsAddRead(&rssOut[2],newl[1],readQual[1],mee[1])) ErrorMsgExit("failed adding read to read set!");
			continue;
		}

    }
    
    kseq_destroy(seq1);
    gzclose(fp[0]);
    gzclose(fpout[0]);
    if (!readSetStatsUpdate(&rssIn[0])) ErrorMsgExit("failed readset update!");
    if (!readSetStatsUpdate(&rssOut[0])) ErrorMsgExit("failed readset update!");    
    if (nInputFastqs==2) {
		kseq_destroy(seq2);
        gzclose(fp[1]);
        gzclose(fpout[1]);
		if (!readSetStatsUpdate(&rssIn[1])) ErrorMsgExit("failed readset update!");
		if (!readSetStatsUpdate(&rssOut[1])) ErrorMsgExit("failed readset update!");    
        if (sflag) 
        {
			gzclose(sfp);
			if (!readSetStatsUpdate(&rssOut[2])) ErrorMsgExit("failed readset update!");    
        }
    }
    
    DebugMsg("Done closing all files");
    
    readSetStatsPrintHeader();     
    readSetStatsPrint(&rssIn[0],&rssIn[0],inputFastqs[0]);
    readSetStatsPrint(&rssOut[0],&rssIn[0],outputFastqs[0]);
    if (nInputFastqs==2)
    {
	readSetStatsPrint(&rssIn[1],&rssIn[1],inputFastqs[0]);
	readSetStatsPrint(&rssOut[1],&rssIn[1],outputFastqs[1]);
	if (sflag) readSetStatsPrint(&rssOut[2],&rssIn[0],singleEndFastqOutFile);

    }
    
    return 1;
}
