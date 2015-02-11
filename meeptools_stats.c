#include "meeptools_stats.h"


static int meeptools_stats_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   meeptools stats [options] \n\n");
    fprintf(stderr, "Options: -f FILE  FASTQ FILE\n");
    fprintf(stderr, "         -?       help\n");
    fprintf(stderr, "\n");
    return 1;
}

int meeptools_stats(int argc, char *argv[])
{
    int i,c,l,l70q20,l70q25,l70q30;
    int fflag=0;
    int debug=0;
    char *fastqFilename;
    gzFile fp;
    kseq_t *seq;
    double readQual;
    double mee;
    double meep;
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
        
    while ((c = getopt(argc, argv, "df:")) != -1)
    {
        switch (c)
        {
        case 'd':
            debug = 1;
            break;
        case 'f':
            fflag = 1;
            fastqFilename = strdup(optarg);
            break;
        case '?':
            return meeptools_stats_usage();
        default:
            return meeptools_stats_usage();
        }
    }
    if (fflag==0){
        ErrorMsg("missing -f option");
        return meeptools_stats_usage();
    }
    if (file_exists(fastqFilename)!=1) {
        sprintf(msg, "%s file does not exist!",fastqFilename);
        ErrorMsgExit(msg);
    }
    
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
            fprintf(stderr,"[%s]: Invalid sequence detected!\n",__func__);
            seq_is_invalid(seq,fastqFilename);
        }

        mee = calculate_mee(seq);

        meep = mee*100.0/l;

        readQual = calculate_qscore_extra(seq,&l70q20,&l70q25,&l70q30);
        
        if (!readSetStatsAddRead(&allReadSetStats[RSALL],l,readQual,mee)) 
        {
            fprintf(stderr,"[%s]: failed adding read to read set %d!\n",__func__,RSALL);
            exit(1);
        }
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
    
    
    fprintf(stdout,"\nNreads\tPercent_reads\tNbases\tPercent_bases\tminRL\tmaxRL\tavgRL\tavgRQ\toverallMEEP\tDescription\n\n");
    for (i=0;i<RSTOTAL;i++)
    {
        fprintf(stdout,"%llu\t%.2f\t%llu\t%.2f\t%u\t%u\t%.2f\t%.2f\t%.4f\t%s\n", \
        allReadSetStats[i].nreads, \
        (allReadSetStats[i].nreads*1.0/allReadSetStats[RSALL].nreads)*100.0, \
        allReadSetStats[i].nbases, \
        (allReadSetStats[i].nbases*1.0/allReadSetStats[RSALL].nbases)*100.0, \
        allReadSetStats[i].minRL, \
        allReadSetStats[i].maxRL, \
        allReadSetStats[i].avgRL, \
        allReadSetStats[i].avgRQ, \
        allReadSetStats[i].overallMEEP, \
        readSetDescriptions[i]);
    }    
    
    
    return 1;
}


