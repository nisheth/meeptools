#include "meeptools_sort.h"
static int meeptools_sort_usage(int extra)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   meeptools sort [options] \n\n Sorts reads by increasing MEEP scores (best to worse)\n\n");
    fprintf(stderr, "Options: -i FILE         FASTQ FILE(S) ... for paired-end read1 and read2 files should be comma separated\n");
    fprintf(stderr, "         -o FILE         FASTQ FILE(S) ... for paired-end read1 and read2 files should be comma separated\n");
    fprintf(stderr, "Optional:\n");
    fprintf(stderr, "         -c FLOAT        MEEP score cut off (between 0 and 20)\n");
    fprintf(stderr, "         -f INTEGER      Offset 33(default) or 64\n");
    fprintf(stderr, "         -t INTEGER      Truncate number of reads in output files\n");
    fprintf(stderr, "         -l INTEGER      read length cut off (default 35)\n");
    fprintf(stderr, "         -q              Append read quality to read comment\n");
    fprintf(stderr, "         -m              Append MEE to read comment\n");
    fprintf(stderr, "         -h              help\n");
    fprintf(stderr, "\n");
    if (extra)
    {
        fprintf(stderr,"Examples:\n");
        fprintf(stderr,"meeptools sort -i a.fastq.gz -o a_meepsorted.fastq.gz -m 1.0\n");
        fprintf(stderr,"meeptools sort -i read1.fastq,read2.fastq -o read1_meep_sorted.fastq.gz,read2_meep_sorted.fastq.gz -m 1.0 -f 64\n");
        fprintf(stderr,"meeptools sort -t 1000 -i read1.fastq,read2.fastq -o read1_meep_sorted.fastq.gz,read2_meep_sorted.fastq.gz\n");
    }
    return 1;
}
int meeptools_sort(int argc, char *argv[])
{
    int z,i,offset=33;
    int l[2];
    int iflag=0;
    int tflag=0;
    int nreadscut=0;
    int oflag=0;
    int cflag=0;
    int mflag=0;
    int qflag=0;
    int lcut=DEFAULTLCUT;
    double mcut=0.0;
    char *inputFastqstr;
    char **inputFastqs;
    int nInputFastqs;
    char *outputFastqstr;
    char **outputFastqs;
    int nOutputFastqs;
    gzFile fp[2],fpout[2];
    kseq_t *seq1,*seq2;
    double mee[2];
    double meep[2];
    double readQual[2];
    readSetStats rssIn[2];
    readSetStats rssOut[2];
    int nreads;
    mykseq *ks,*allReads=NULL;    
    char msg[STDSTRLEN];
    
    for (i=0;i<2;i++) 
    {
        mee[i]=-1.0;
        meep[i]=-1.0;
        readQual[i]=0.0;
        fpout[i]=NULL;
        if (!readSetStatsInit(&rssIn[i])) ErrorMsgExit("failed initialization of read set!");
        if (!readSetStatsInit(&rssOut[i])) ErrorMsgExit("failed initialization of read set!");
    }
    
    while ((z = getopt(argc, argv, "i:f:o:c:l:t:mqh")) != -1)
    {
        switch (z)
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
        case 't':
            tflag = 1;
            nreadscut = atoi(optarg);
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
            return meeptools_sort_usage(1);
        default:
            return meeptools_sort_usage(0);
        }
    }

    sprintf(msg,"MEEP cut off = %.2f",mcut);
    DebugMsg(msg);
    sprintf(msg,"RL cut off = %d",lcut);
    DebugMsg(msg);
    if (tflag) {
    sprintf(msg,"Truncate after %d reads",nreadscut);
    DebugMsg(msg);
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
        ErrorMsg("missing -f option");
        return meeptools_sort_usage(0);
    }
    if (oflag==0)
    {
        ErrorMsg("missing -o option");
        return meeptools_sort_usage(0);
    }
    if (lcut<0)
    {
        ErrorMsg("Invalid -l option");
        return meeptools_sort_usage(0);
    }
    if ((cflag==1) && (mcut<0 || mcut>20))
    {
        ErrorMsg("MEEP cut off value out of bounds!");
        return meeptools_sort_usage(0);
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
        
    }
    
    while (0==0)
    {
        l[0] = kseq_read(seq1);
        if (l[0] == -1) break;
        if (l[0] == -2)
        {
            ErrorMsg("Invalid sequence detected!");
            seqIsInvalid(seq1,inputFastqs[0]);
        }
		
        mee[0] = seqCalculateMEE(seq1);
        meep[0] = mee[0]*100.0/l[0];
        readQual[0] = seqCalculateQScore(seq1,offset);
	
        if (!readSetStatsAddRead(&rssIn[0],l[0],readQual[0],mee[0])) ErrorMsgExit("failed adding read to read set!");
	
		if (nInputFastqs==2)
        {
            l[1] = kseq_read(seq2);
            if (l[1] == -1) break;
            if (l[1] == -2)
            {
                ErrorMsg("Invalid sequence detected!");
                seqIsInvalid(seq2,inputFastqs[1]);
            }            
            
            mee[1] = seqCalculateMEE(seq2);
            meep[1] = mee[1]*100.0/l[1];
	    	readQual[1] = seqCalculateQScore(seq2,offset);
	    	
	    	if (!readSetStatsAddRead(&rssIn[1],l[1],readQual[1],mee[1])) ErrorMsgExit("failed adding read to read set!");
        }
	
    	if(nreads%PRINTINTERVAL == 0)
		{
			sprintf(msg,"%d reads processed.",nreads);
			PrintMsg(msg);
		}

        nreads++;
        ks=(mykseq*)malloc(sizeof(mykseq));
        ks->name1=(char *)malloc(sizeof(char)*(seq1->name.l+1));
        memset(ks->name1, '\0', sizeof(seq1->name.l+1));
        strcpy(ks->name1, seq1->name.s);
        ks->comment1=(char *)malloc(sizeof(char)*(seq1->comment.l+1));
        memset(ks->comment1, '\0', sizeof(seq1->comment.l+1));
        strcpy(ks->comment1, seq1->comment.s);
        ks->bases1=(char *)malloc(sizeof(char)*(seq1->seq.l+1));
        memset(ks->bases1, '\0', sizeof(seq1->seq.l+1));
        strcpy(ks->bases1, seq1->seq.s);
        ks->quality1=(char *)malloc(sizeof(char)*(seq1->qual.l+1));
        memset(ks->quality1, '\0', sizeof(seq1->qual.l+1));
        strcpy(ks->quality1, seq1->qual.s);
        ks->mee1=mee[0];
        ks->l1=l[0];
        ks->meep1=meep[0];
        ks->rq1=readQual[0];
        ks->combinedMEEP=meep[0];
		
		if (nInputFastqs==2) 
		{
		    ks->name2=(char *)malloc(sizeof(char)*(seq2->name.l+1));
			memset(ks->name2, '\0', sizeof(seq2->name.l+1));
			strcpy(ks->name2, seq2->name.s);
			ks->comment2=(char *)malloc(sizeof(char)*(seq2->comment.l+1));
			memset(ks->comment2, '\0', sizeof(seq2->comment.l+1));
			strcpy(ks->comment2, seq2->comment.s);
			ks->bases2=(char *)malloc(sizeof(char)*(seq2->seq.l+1));
			memset(ks->bases2, '\0', sizeof(seq2->seq.l+1));
			strcpy(ks->bases2, seq2->seq.s);
			ks->quality2=(char *)malloc(sizeof(char)*(seq2->qual.l+1));
			memset(ks->quality2, '\0', sizeof(seq2->qual.l+1));
			strcpy(ks->quality2, seq2->qual.s);
			ks->mee2=mee[1];
			ks->l2=l[1];
			ks->meep2=meep[1];
			ks->rq2=readQual[1];
			ks->combinedMEEP=(mee[0]+mee[1])*100.0/(l[0]+l[1]);
		}
		
		HASH_ADD_INT(allReads,count,ks);
	}
	
	kseq_destroy(seq1);
    gzclose(fp[0]);
    if (!readSetStatsUpdate(&rssIn[0])) ErrorMsgExit("failed readset update!");
    if (nInputFastqs==2) {
		kseq_destroy(seq2);
        gzclose(fp[1]);
		if (!readSetStatsUpdate(&rssIn[1])) ErrorMsgExit("failed readset update!"); 
    }
    
    DebugMsg("Done closing all input files");
    
    DebugMsg("Sorting now...");
    HASH_SORT(allReads, mykseqMEEPSort);
	DebugMsg("Done sorting!");
    	
	i=0;
	for(ks=allReads; ks != NULL; ks=ks->hh.next) {
		if (cflag && (ks->combinedMEEP>mcut)) break;
		if (ks->l1<lcut) continue;
		if ((nInputFastqs==2) && (ks->l2<lcut)) continue;
		i++;
		if (!readSetStatsAddRead(&rssOut[0],ks->l1,ks->rq1,ks->mee1)) ErrorMsgExit("failed adding read to read set!");
		if (nInputFastqs==2) {
			if (!readSetStatsAddRead(&rssOut[1],ks->l2,ks->rq2,ks->mee2)) ErrorMsgExit("failed adding read to read set!");
		}
		if(!mykseqWriteToFile(ks,fpout,mflag,qflag)) ErrorMsgExit("Failed to write to output files");
		if (tflag && i==nreadscut) break;

	}

    gzclose(fpout[0]);
    if (!readSetStatsUpdate(&rssOut[0])) ErrorMsgExit("failed readset update!");    
    if (nInputFastqs==2) {
        gzclose(fpout[1]);
		if (!readSetStatsUpdate(&rssOut[1])) ErrorMsgExit("failed readset update!");    
    }
    
    DebugMsg("Done closing all output files");

    readSetStatsPrintHeader();     
    readSetStatsPrint(&rssIn[0],&rssIn[0],inputFastqs[0]);
    readSetStatsPrint(&rssOut[0],&rssIn[0],outputFastqs[0]);
    if (nInputFastqs==2)
    {
	readSetStatsPrint(&rssIn[1],&rssIn[1],inputFastqs[0]);
	readSetStatsPrint(&rssOut[1],&rssIn[1],outputFastqs[1]);
    }
    
    
	return 1;
        
}

int mykseqMEEPSort(mykseq *a, mykseq *b) {
	if (a->combinedMEEP < b->combinedMEEP) {
		return -1;
	} else if (a->combinedMEEP > b->combinedMEEP) {
		return 1;
	} else {
		return 0;
	}
}

