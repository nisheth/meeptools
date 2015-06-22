#include "meeptools_stats.h"


static int meeptools_stats_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   meeptools stats [options] \n\n Generates MEEP score and other stats for fastq file.\n\n");
    fprintf(stderr, "Options: -i FILE  FASTQ FILE\n");
    fprintf(stderr, "         -s       Output subset stats also.\n");
    fprintf(stderr, "         -f INT   Offset 33(default) or 64 (Optional)\n.");
    fprintf(stderr, "         -h       help\n");
    fprintf(stderr, "\n");
    return 1;
}

int meeptools_stats(int argc, char *argv[])
{
    int i,c,l,l70q20,l70q25,l70q30,offset=33;
    int iflag=0;
    int sflag=0;
    int jflag=0;
    int dflag=0;
    char *fastqFilename;
    char jsonFilename[500];
    FILE *jf;
    gzFile fp;
    kseq_t *seq;
    double readQual=0.0;
    double mee=0.0;
    double meep=0.0;
    char msg[STDSTRLEN];
    
    char* readSetDescriptions[] = {RS_DESC};
        
    readSetStats allReadSetStats[RSTOTAL];
    
    
    for (i=0;i<RSTOTAL;i++){
        if (!readSetStatsInit(&allReadSetStats[i]))
        {
            sprintf(msg,"failed initialization of read set %s!",readSetDescriptions[i]);
            ErrorMsgExit(msg);
        }
    }
        
    while ((c = getopt(argc, argv, "i:f:sjdh?")) != -1)
    {
        switch (c)
        {
        case 'i':
            iflag = 1;
            fastqFilename = strdup(optarg);
            break;
        case 's':
            sflag = 1;
            break;
        case 'd':
            dflag = 1;
            break;
        case 'j':
            jflag = 1;
            break;
        case 'f':
            offset = atoi(optarg);
            break;
        case 'h':
            return meeptools_stats_usage();
        default:
            return meeptools_stats_usage();
        }
    }
    if (iflag==0){
        ErrorMsg("missing -i option");
        return meeptools_stats_usage();
    }
    if (file_exists(fastqFilename)!=1) {
        sprintf(msg, "%s file does not exist!",fastqFilename);
        ErrorMsgExit(msg);
    }
    if (offset !=33 && offset !=64 ) {
    	ErrorMsgExit("Offset value can only be 33(default) or 64!");
    }
    if (!init_q2mee_hash(offset))
    {
        ErrorMsgExit("Initialization of q2mee failed!");
    }
	DebugMsg("Initializing q2mee hash complete.");
    
    if ((fp = gzopen(fastqFilename, "rb")) == NULL )
    {
        sprintf(msg, "%s file cannot be opened!",fastqFilename);
        ErrorMsgExit(msg);    
    } 
    seq = kseq_init(fp);

    sprintf(msg,"FASTQ file %s opened.",fastqFilename);
    DebugMsg(msg);

    while (0==0)
    {
        l = kseq_read(seq);
        if (l == -1) break;
        if (l == -2) {
            ErrorMsg("Invalid sequence detected!");
            seqIsInvalid(seq,fastqFilename);
        }

        mee = seqCalculateMEE(seq);

        meep = mee*100.0/l;

        readQual = seqCalculateQScoreExtra(seq,offset,&l70q20,&l70q25,&l70q30);
        
        if (!readSetStatsAddRead(&allReadSetStats[RSALL],l,readQual,mee)) 
        {
            fprintf(stderr,"[%s]: failed adding read to read set %d!\n",__func__,RSALL);
            exit(1);
        }
        if (sflag)
        {
            if (l70q20)
            {
                if (!readSetStatsAddRead(&allReadSetStats[RSL70Q20],l,readQual,mee)) 
                {
                    fprintf(stderr,"[%s]: failed adding read to read set %d!\n",__func__,RSL70Q20);
                    exit(1);
                }    
            }
            if (l70q25)
            {
                if (!readSetStatsAddRead(&allReadSetStats[RSL70Q25],l,readQual,mee)) 
                {
                    fprintf(stderr,"[%s]: failed adding read to read set %d!\n",__func__,RSL70Q25);
                    exit(1);
                }    
            }
            if (l70q30)
            {
                if (!readSetStatsAddRead(&allReadSetStats[RSL70Q30],l,readQual,mee)) 
                {
                    fprintf(stderr,"[%s]: failed adding read to read set %d!\n",__func__,RSL70Q30);
                    exit(1);
                }    
            }
            if (readQual>=20)
            {
                if (!readSetStatsAddRead(&allReadSetStats[RSQ20],l,readQual,mee)) 
                {
                    fprintf(stderr,"[%s]: failed adding read to read set %d!\n",__func__,RSQ20);
                    exit(1);
                }    
            }
            if (readQual>=25)
            {
                if (!readSetStatsAddRead(&allReadSetStats[RSQ25],l,readQual,mee)) 
                {
                    fprintf(stderr,"[%s]: failed adding read to read set %d!\n",__func__,RSQ25);
                    exit(1);
                }    
            }
            if (readQual>=30)
            {
                if (!readSetStatsAddRead(&allReadSetStats[RSQ30],l,readQual,mee)) 
                {
                    fprintf(stderr,"[%s]: failed adding read to read set %d!\n",__func__,RSQ30);
                    exit(1);
                }    
            }
            if (meep<1)
            {
                if (!readSetStatsAddRead(&allReadSetStats[RSMEEP1],l,readQual,mee)) 
                {
                    fprintf(stderr,"[%s]: failed adding read to read set %d!\n",__func__,RSMEEP1);
                    exit(1);
                }    
            }
            if (meep<2)
            {
                if (!readSetStatsAddRead(&allReadSetStats[RSMEEP2],l,readQual,mee)) 
                {
                    fprintf(stderr,"[%s]: failed adding read to read set %d!\n",__func__,RSMEEP2);
                    exit(1);
                }    
            }
        }
    }
    
    kseq_destroy(seq);
    gzclose(fp);
    
    sprintf(msg,"FASTQ file %s closed.",fastqFilename);
    DebugMsg(msg);
    
    for (i=0;i<RSTOTAL;i++){
        if (!readSetStatsUpdate(&allReadSetStats[i]))
        {
            fprintf(stderr,"[%s]: update failed of read set %d!\n",__func__,i);
            exit(1);
        }
    }
    

    readSetStatsPrintHeader();    
    if (sflag)
    {
        c=RSTOTAL;
    }
    else
    {
        c=1;
    }
    for (i=0;i<c;i++)
    {
        readSetStatsPrint(&allReadSetStats[i],&allReadSetStats[RSALL],readSetDescriptions[i]);
    }    
    
    if (dflag) {
        readSetStatsPrintHist(&allReadSetStats[RSALL],fastqFilename);
	}
	
	if (jflag) {
		//sprintf(jsonFilename,"%s.json",fastqFilename);
		//jf=fopen(jsonFilename,"w");
        readSetStatsPrintJson(&allReadSetStats[RSALL],fastqFilename);
        //fclose(jf);
    }
    return 1;
}



