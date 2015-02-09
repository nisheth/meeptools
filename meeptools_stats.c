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
    int c;
    int fflag=0;
    int debug=0;
    char *fastqFilename;
    
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
        fprintf(stderr, "[meeptools_stats]: missing -f option\n");
        return meeptools_stats_usage();
    }
    if (file_exists(fastqFilename)!=1) {
        fprintf(stderr, "[meeptools_stats]: %s file does not exist!\n",fastqFilename);
        exit(1);
    } else {
        fprintf(stdout, "[meeptools_stats]: %s file found\n",fastqFilename);
        exit(0);
    }
}
