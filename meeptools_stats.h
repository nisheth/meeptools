#include "utils.h"

#define RSALL 0
#define RSQ20 1
#define RSQ25 2
#define RSQ30 3
#define RSL70Q20 4
#define RSL70Q25 5
#define RSL70Q30 6
#define RSMEEP1 7
#define RSMEEP2 8
#define RSTOTAL 9

#define RS_DESC "All_reads", \
"Reads_with_avgRQ_atleast_Q20", \
"Reads_with_avgRQ_atleast_Q25", \
"Reads_with_avgRQ_atleast_Q30", \
"Reads_with_70%_bases_atleast_Q20", \
"Reads_with_70%_bases_atleast_Q25", \
"Reads_with_70%_bases_atleast_Q30", \
"Reads_with_MEEP_LT_1", \
"Reads_with_MEEP_LT_2"


static int meeptools_stats_usage();

