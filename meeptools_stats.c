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
    int c,l,l70q20,l70q25,l70q30;
    unsigned long long int nreadsl70q20=0;
    unsigned long long int nreadsl70q25=0;
    unsigned long long int nreadsl70q30=0;
    int fflag=0;
    int debug=0;
    char *fastqFilename;
    gzFile fp;
    kseq_t *seq;
    unsigned long long int nreads=0;
    unsigned long long int nbases=0;
    unsigned long long int nreadsQ20=0;
    unsigned long long int nreadsQ25=0;
    unsigned long long int nreadsQ30=0;
    unsigned long long int meep1=0;
    unsigned long long int meep2=0;
    unsigned long long int meep1Nbases=0;
    unsigned long long int meep2Nbases=0;
    double avgReadLen=0.0;
    double avgReadQual=0.0;
    double avgReadQualMeep1=0.0;
    double avgReadQualMeep2=0.0;
    double readQual;
    double readQualSum=0.0;
    double readQualSumMeep1=0.0;
    double readQualSumMeep2=0.0;
    double mee;
    double meep;
    double overallMeep=0;
    double overallMeep1=0;
    double overallMeep2=0;
    
    unsigned int minRL=(unsigned int)MAXREADLENGTH;
    unsigned int maxRL=0;
    double preads_meep1=0.0;
    double preads_meep2=0.0;
    double preads_Q20=0.0;
    double preads_Q25=0.0;
    double preads_Q30=0.0;
    double preads_l70q20=0.0;
    double preads_l70q25=0.0;
    double preads_l70q30=0.0;
    
    
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
        fprintf(stderr, "[%s]: missing -f option\n",__func__);
        return meeptools_stats_usage();
    }
    if (file_exists(fastqFilename)!=1) {
        fprintf(stderr, "[%s]: %s file does not exist!\n",__func__,fastqFilename);
        exit(1);
    }
    
    fp = gzopen(fastqFilename, "rb"); 
    seq = kseq_init(fp);

    while (0==0)
    {
        l = kseq_read(seq);
        if (l == -1) break;
        if (l == -2) {
            fprintf(stderr,"[%s]: Invalid sequence detected!\n",__func__);
            seq_is_invalid(seq,fastqFilename);
        }
        nreads++;
        nbases+=l;
        if (l>maxRL) maxRL=l;
        if (l<minRL) minRL=l;

        mee = calculate_mee(seq);

        meep = mee*100.0/l;
        overallMeep += mee;

        readQual = calculate_qscore_extra(seq,&l70q20,&l70q25,&l70q30);

        readQualSum += readQual;
        nreadsl70q20 += l70q20;
        nreadsl70q25 += l70q25;
        nreadsl70q30 += l70q30;
        
        if(readQual>=20) nreadsQ20++;
        if(readQual>=25) nreadsQ25++;
        if(readQual>=30) nreadsQ30++;        
        
        if (meep<2) {
            meep2++;
            meep2Nbases += l;
            overallMeep2 += mee;
            readQualSumMeep2 += readQual;
            if (meep<1) {
                meep1++;
                meep1Nbases += l;
                overallMeep1 += mee;
                readQualSumMeep1 += readQual;
            }
        }


    }
    kseq_destroy(seq);
    gzclose(fp);
    if (nreads>0) {
        avgReadQual=readQualSum/nreads;
        avgReadLen=nbases*1.0/nreads;
        preads_meep1=meep1*100.0/nreads;
        preads_meep2=meep2*100.0/nreads;
        overallMeep=overallMeep*100.0/nbases;
        preads_l70q20=nreadsl70q20*100.0/nreads;
        preads_l70q25=nreadsl70q25*100.0/nreads;
        preads_l70q30=nreadsl70q30*100.0/nreads;
        preads_Q20=nreadsQ20*100/nreads;
        preads_Q25=nreadsQ25*100/nreads;
        preads_Q30=nreadsQ30*100/nreads;
    }
    if (meep1Nbases>0) {
        overallMeep1=overallMeep1*100.0/meep1Nbases;
        avgReadQualMeep1=readQualSumMeep1/meep1;
    }
    if (meep2Nbases>0) {
        overallMeep2=overallMeep2*100.0/meep2Nbases;
        avgReadQualMeep2=readQualSumMeep2/meep2;
    }
    
    printf("Number_of_reads=%llu\n",nreads);
    printf("Number_of_bases=%llu\n",nbases);
    printf("Minimum_read_length=%d\n",minRL);
    printf("Maximum_read_length=%d\n",maxRL);
    printf("Average_read_length=%.2f\n",avgReadLen);
    printf("Average_read_quality=%.2f\n",avgReadQual);
    printf("Number_of_reads_with_avgReadQuality_GT_Q20=%llu\n",nreadsQ20);
    printf("Percent_of_reads_with_avgReadQuality_GT_Q20=%.2f\n",preads_Q20);
    printf("Number_of_reads_with_avgReadQuality_GT_Q25=%llu\n",nreadsQ25);
    printf("Percent_of_reads_with_avgReadQuality_GT_Q25=%.2f\n",preads_Q25);
    printf("Number_of_reads_with_avgReadQuality_GT_Q30=%llu\n",nreadsQ30);
    printf("Percent_of_reads_with_avgReadQuality_GT_Q30=%.2f\n",preads_Q30);
    printf("Number_of_reads_with_70_percent_of_bases_GT_Q20=%llu\n",nreadsl70q20);
    printf("Percent_of_reads_with_70_percent_of_bases_GT_Q20=%.2f\n",preads_l70q20);
    printf("Number_of_reads_with_70_percent_of_bases_GT_Q25=%llu\n",nreadsl70q25);
    printf("Percent_of_reads_with_70_percent_of_bases_GT_Q25=%.2f\n",preads_l70q25);
    printf("Number_of_reads_with_70_percent_of_bases_GT_Q30=%llu\n",nreadsl70q30);
    printf("Percent_of_reads_with_70_percent_of_bases_GT_Q30=%.2f\n",preads_l70q30);
    printf("Overall_MEEP=%.4f\n",overallMeep);
    printf("Number_of_MEEP1_reads=%llu\n",meep1);
    printf("Percent_MEEP1_reads=%.2f\n",preads_meep1);
    printf("Overall_MEEP_of_MEEP1_reads=%.4f\n",overallMeep1);
    printf("Average_read_quality_of_MEEP1_reads=%.2f\n",avgReadQualMeep1);
    printf("Number_of_MEEP2_reads=%llu\n",meep2);
    printf("Percent_MEEP2_reads=%.2f\n",preads_meep2);
    printf("Overall_MEEP_of_MEEP2_reads=%.4f\n",overallMeep2);
    printf("Average_read_quality_of_MEEP2_reads=%.2f\n",avgReadQualMeep2);
    
    return 1;
}


