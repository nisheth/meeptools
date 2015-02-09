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
#define COMMONBLOCK "1"
#endif

int meeptools_append(int argc, char *argv[]);
int meeptools_filter(int argc, char *argv[]);
int meeptools_sort(int argc, char *argv[]);
int meeptools_stats(int argc, char *argv[]);
int meeptools_subset(int argc, char *argv[]);
int meeptools_trim(int argc, char *argv[]);
