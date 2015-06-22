#include "meeptools.h"

static int usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: meeptools (Tools to calculate maximum expected error in a FASTQ read as a percentage of read length)\n");
    fprintf(stderr, "Version: %s\n\n", PACKAGE_VERSION);
    fprintf(stderr, "Usage:   meeptools <command> [options]\n\n");
    fprintf(stderr, "Command: append      append MEEP score to FASTQ file\n");
    fprintf(stderr, "         filter      filter FASTQ file based on MEEP score\n");
    fprintf(stderr, "         sort        sort reads by MEEP score\n");
    fprintf(stderr, "         stats       MEEP score based stats for FASTQ file\n");
    fprintf(stderr, "         trim        trim reads based on MEEP score\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[])
{

    if (argc < 2) return usage();
    
    if (strcmp(argv[1], "append") == 0) return meeptools_append(argc-1, argv+1);
    else if (strcmp(argv[1], "filter") == 0) return meeptools_filter(argc-1, argv+1);
    else if (strcmp(argv[1], "sort") == 0) return meeptools_sort(argc-1, argv+1);
    else if (strcmp(argv[1], "stats") == 0) return meeptools_stats(argc-1, argv+1);
    else if (strcmp(argv[1], "trim") == 0) return meeptools_trim(argc-1, argv+1);
    else {
            fprintf(stderr, "[%s] unrecognized command '%s'\n", __func__,argv[1]);
            return 1;
    }

    return 1;

}
