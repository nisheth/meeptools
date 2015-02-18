#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "uthash.h"
#include "kseq.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1"
#endif

#ifndef COMMONBLOCK
#define COMMONBLOCK 1
#endif

#ifndef MAXREADLENGTH
#define MAXREADLENGTH 10000
#endif

#ifndef STDSTRLEN
#define STDSTRLEN 256
#endif

#ifndef PRINTINTERVAL
#define PRINTINTERVAL 10000
#endif

#ifndef DEFAULTLCUT
#define DEFAULTLCUT 35
#endif

#ifdef DEBUG
#define DebugMsg(msg) fprintf(stderr,"[%s]: %s\n",__func__,msg)
#else
#define DebugMsg(msg) NULL
#endif

#define ErrorMsg(msg) fprintf(stderr,"ERROR [%s]: %s\n",__func__,msg)
#define PrintMsg(msg) fprintf(stderr,"MESSAGE [%s]: %s\n",__func__,msg)
#define ErrorMsgExit(msg) {fprintf(stderr,"ERROR [%s]: %s\n",__func__,msg);exit(1);}

struct meeHash{
    int Q;//key
    double MEE;
    UT_hash_handle hh;
};


int meeptools_append(int argc, char *argv[]);
int meeptools_filter(int argc, char *argv[]);
int meeptools_sort(int argc, char *argv[]);
int meeptools_stats(int argc, char *argv[]);
int meeptools_subset(int argc, char *argv[]);
int meeptools_trim(int argc, char *argv[]);

KSEQ_INIT(gzFile, gzread)
