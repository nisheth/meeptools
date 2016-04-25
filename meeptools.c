#include "meeptools.h"

void bseq1_calculate_MEE_Qsum(bseq1_t *seq,double *q2mee)
{
    int i;
    int j;
    seq->o_mee=0.0;
    seq->o_qsum=0;
    seq->o_meep=0.0;
    seq->o_avgq=0.0;
    for (i=0;i<seq->l_seq;i++) {
    	j=(int)seq->qual[i];
    	seq->o_qsum+=j;
        seq->o_mee+=q2mee[j];
    }
    seq->o_meep=seq->o_mee*100.0/seq->l_seq;
    seq->o_qsum-=(33 * seq->l_seq);
    seq->o_avgq=(seq->o_qsum*1.0)/seq->l_seq;
// 					fprintf(stderr,"[%s]%d\t%d\t%.4f\t%.4f\n",__func__,seq->trim_start,seq->trim_len,seq->mee,seq->meep);
}

void bseq1_copy_MEE_Qsum(bseq1_t *seq)
{
    seq->mee=seq->o_mee;
    seq->qsum=seq->o_qsum;
    seq->meep=seq->o_meep;
    seq->avgq=seq->o_avgq;
    seq->trim_start=0;
    seq->trim_len=seq->l_seq;
// 					fprintf(stderr,"[%s]%d\t%d\t%.4f\t%.4f\n",__func__,seq->trim_start,seq->trim_len,seq->mee,seq->meep);
}


void bseq1_update_MEE_Qsum(bseq1_t *seq,double *q2mee)
{
    int i;
    int j;
    seq->mee=0.0;
    seq->qsum=0;
    seq->meep=0.0;
    seq->avgq=0.0;
    for (i=seq->trim_start;i<(seq->trim_start+seq->trim_len);i++) {
    	j=(int)seq->qual[i];
    	seq->qsum+=j;
        seq->mee+=q2mee[j];
    }
    seq->meep=seq->mee*100.0/seq->trim_len;
    seq->qsum-=(33 * seq->trim_len);
    seq->avgq=(seq->qsum*1.0)/seq->trim_len;
// 					fprintf(stderr,"[%s]%d\t%d\t%.4f\t%.4f\n",__func__,seq->trim_start,seq->trim_len,seq->mee,seq->meep);
}

int bseq1_MEE_trim_fast(bseq1_t *seq,double *q2mee,int lcut,double mcut)
{
    int i;
    int j;
    int k;
    int l2;
    double x,y;
    int l=seq->trim_len - seq->trim_start;
    double mee[l];
    double sumatlen[l];
    double previous;
    
    if (seq->trim_len < lcut) return 0;
    if (seq->meep < mcut) return 1;

    for (i=seq->trim_start;i<seq->trim_len;i++) 
    {
    	j=(int)seq->qual[i];
    	mee[i]=q2mee[j];
    }
	for (i=0;i<=(l-lcut);i++) {
		for (j=l;j>=lcut;j--){
			if (j-i<lcut) break;
			if (j==l) {
				x=0.0;
				for (k=i;k<j;k++) x+=mee[k];
				previous = x;
			} else {
				previous -= mee[j+1];
				x = previous;
			}
			if ((x*100)/(j-i) < mcut) {
				seq->trim_len=j-i;
				seq->trim_start=i;
				bseq1_update_MEE_Qsum(seq,q2mee);
				return 1;
			}
		}
	}
	
	return 0; // no amount of trimming fit the mcut threshold
}

int bseq1_MEE_trim_slow(bseq1_t *seq,double *q2mee,int lcut,double mcut)
{
    int i;
    int j;
    int k;
    int l2;
    double x;
    int l=seq->trim_len - seq->trim_start;
    double mee[l];
    double sumatlen[l];
    double previous;
    
    if (seq->l_seq < lcut) return 0;
    if (seq->meep < mcut) return 1;

    for (i=seq->trim_start;i<seq->trim_len;i++) 
    {
    	j=(int)seq->qual[i];
    	mee[i]=q2mee[j];
    }
	

	for (j=l;j>=lcut;j--){
	for (i=0;i<=(l-lcut);i++) {
			if (i+j>l) break;
				x=0.0;
				for (k=i;k<i+j;k++) x+=mee[k];
			if ((x*100)/j < mcut) {
				seq->trim_len=j;
				seq->trim_start=i;
				bseq1_update_MEE_Qsum(seq,q2mee);
				return 1;
			}
		}
	}	

	return 0; // no amount of trimming fit the mcut threshold
}

int bseq1_MEE_trim(bseq1_t *seq,double *q2mee,int lcut,double mcut)
{
    int i;
    int j;
    int k;
    int l2;
    double x;
    int l=seq->trim_len - seq->trim_start;
    double mee[l];
    double sumatlen[l];
    double previous;
    
    if (seq->l_seq < lcut) return 0;
    if (seq->meep < mcut) return 1;

    for (i=seq->trim_start;i<seq->trim_len;i++) 
    {
    	j=(int)seq->qual[i];
    	mee[i]=q2mee[j];
    }
	

	for (i=0;i<=(l-lcut);i++) {
		for (j=l;j>=lcut;j--){
			if (j-i<lcut) break;
				x=0.0;
				for (k=i;k<j;k++) x+=mee[k];
			if ((x*100)/(j-i) < mcut) {
				seq->trim_len=j-i;
				seq->trim_start=i;
				bseq1_update_MEE_Qsum(seq,q2mee);
				return 1;
			}
		}
	}
	return 0; // no amount of trimming fit the mcut threshold
}
void opt_init(opt_t *opt)
{
	int i;
	memset(opt, 0, sizeof(opt_t));
	opt->chunk_size = 10000; // number of reads read in one chunk
	opt->n_threads = 4; // default number of threads is 4
	for (i=0;i<33;i++) { opt->q2mee[i]=0.0; }
	for (i=33;i<75;i++) {
		opt->q2mee[i]=pow(10.0,(-1.0*(i-33)/10));
	}
}

// for callback function
// This used to be call for each read, now called for each chunk of reads

static void kt_for_cb(void *_data, long k, int tid) 
{
	//fprintf(stderr,"[%s] HERE tid=%d\n",__func__,tid);
	step_t *data = (step_t*)_data;
	shared_t *es = data->es;
	bseq1_t *s = (bseq1_t *)data->seqsarray[k];
	bseq1_t *s2 = (bseq1_t *)NULL;
	double *q2mee = (double *)es->opt->q2mee;
	int i,result;

	for (i=0;i<data->n_seqs[k];i++) {
		bseq1_calculate_MEE_Qsum(&s[i],q2mee);
		// fprintf(stderr,"[%s] i=%d mee=%f\n",__func__,i,s[i].o_mee);
		s[i].output_mode=NO_OUTPUT;
		if (es->meeptools_mode==MMODE_APPEND) s[i].output_mode=PAIRED_OUTPUT;
		bseq1_copy_MEE_Qsum(&s[i]);
	}
	if (es->paired==1) {
		s2 = (bseq1_t *)data->seqsarray2[k];
		for (i=0;i<data->n_seqs2[k];i++) {
			bseq1_calculate_MEE_Qsum(&s2[i],q2mee);
			s2[i].output_mode=NO_OUTPUT;
			if (es->meeptools_mode==MMODE_APPEND) s2[i].output_mode=PAIRED_OUTPUT;
			bseq1_copy_MEE_Qsum(&s2[i]);
		}
	}

	if (es->meeptools_mode==MMODE_FILTER) {
		for (i=0;i<data->n_seqs[k];i++) {
			if (s[i].meep > es->mcut) {
				s[i].output_mode = NO_OUTPUT;
			} else {
				s[i].output_mode = PAIRED_OUTPUT; // PAIRED_OUTPUT is misleading term .. just means R1IN ---> R1OUT 
			}
		}
		if (es->paired==1) {
			for (i=0;i<data->n_seqs2[k];i++) {
				if (s2[i].meep > es->mcut) {
					s2[i].output_mode = NO_OUTPUT;
					if (s[i].output_mode == PAIRED_OUTPUT) s[i].output_mode = SINGLE_OUTPUT;
				} else {
					s2[i].output_mode = PAIRED_OUTPUT; // PAIRED_OUTPUT is misleading term .. just means R1IN ---> R1OUT 
					if (s[i].output_mode == NO_OUTPUT) s2[i].output_mode = SINGLE_OUTPUT;
				}
			}
		}
	}

	if (es->meeptools_mode==MMODE_TRIM) {
		for (i=0;i<data->n_seqs[k];i++) {
			if (bseq1_MEE_trim(&s[i],q2mee,es->lcut,es->mcut)==0) {
				s[i].output_mode = NO_OUTPUT;
			} else {
				s[i].output_mode = PAIRED_OUTPUT; // PAIRED_OUTPUT is misleading term .. just means R1IN ---> R1OUT 
			}
		}
		if (es->paired==1) {
			for (i=0;i<data->n_seqs2[k];i++) {
				if (bseq1_MEE_trim(&s2[i],q2mee,es->lcut,es->mcut)==0) {
					s2[i].output_mode = NO_OUTPUT;
					if (s[i].output_mode == PAIRED_OUTPUT) s[i].output_mode = SINGLE_OUTPUT;
				} else {
					s2[i].output_mode = PAIRED_OUTPUT; // PAIRED_OUTPUT is misleading term .. just means R1IN ---> R1OUT 
					if (s[i].output_mode == NO_OUTPUT) s2[i].output_mode = SINGLE_OUTPUT;
				}
			}
		}
	}
	if (es->meeptools_mode==MMODE_TRIM_FAST) {
		for (i=0;i<data->n_seqs[k];i++) {
			if (bseq1_MEE_trim_fast(&s[i],q2mee,es->lcut,es->mcut)==0) {
				s[i].output_mode = NO_OUTPUT;
			} else {
				s[i].output_mode = PAIRED_OUTPUT; // PAIRED_OUTPUT is misleading term .. just means R1IN ---> R1OUT 
			}
		}
		if (es->paired==1) {
			for (i=0;i<data->n_seqs2[k];i++) {
				if (bseq1_MEE_trim_fast(&s2[i],q2mee,es->lcut,es->mcut)==0) {
					s2[i].output_mode = NO_OUTPUT;
					if (s[i].output_mode == PAIRED_OUTPUT) s[i].output_mode = SINGLE_OUTPUT;
				} else {
					s2[i].output_mode = PAIRED_OUTPUT; // PAIRED_OUTPUT is misleading term .. just means R1IN ---> R1OUT 
					if (s[i].output_mode == NO_OUTPUT) s2[i].output_mode = SINGLE_OUTPUT;
				}
			}
		}
	}
	if (es->meeptools_mode==MMODE_TRIM_SLOW) {
		for (i=0;i<data->n_seqs[k];i++) {
			if (bseq1_MEE_trim_slow(&s[i],q2mee,es->lcut,es->mcut)==0) {
				s[i].output_mode = NO_OUTPUT;
			} else {
				s[i].output_mode = PAIRED_OUTPUT; // PAIRED_OUTPUT is misleading term .. just means R1IN ---> R1OUT 
			}
		}
		if (es->paired==1) {
			for (i=0;i<data->n_seqs2[k];i++) {
				if (bseq1_MEE_trim_slow(&s2[i],q2mee,es->lcut,es->mcut)==0) {
					s2[i].output_mode = NO_OUTPUT;
					if (s[i].output_mode == PAIRED_OUTPUT) s[i].output_mode = SINGLE_OUTPUT;
				} else {
					s2[i].output_mode = PAIRED_OUTPUT; // PAIRED_OUTPUT is misleading term .. just means R1IN ---> R1OUT 
					if (s[i].output_mode == NO_OUTPUT) s2[i].output_mode = SINGLE_OUTPUT;
				}
			}
		}
	}
}


// actual pipeline callback function
// step 1, i.e. step=0 ... for reading everything in chunks
// step 2, i.e. step=1 ... for sending different chunks to different cpus and Processing reads
// step 3, i.e. step=2 ... for collecting all processed reads and printing out
void *kt_pipeline_cb(void *shared, int step, void *_data) 
{
	shared_t *es = (shared_t*)shared;
	if (step == 0) { // read
		step_t *ret;
		ret = calloc(1, sizeof(step_t));
		ret->items = 10;
		ret->seqsarray = calloc(ret->items, sizeof(bseq1_t));
		ret->n_seqs = calloc(ret->items,sizeof(int));
		ret->es = es;
		if (es->paired == 1) {
			ret->n_seqs2 = calloc(ret->items,sizeof(int));
			ret->seqsarray2 = calloc(ret->items, sizeof(bseq1_t));
		}
		int keep_comment = 1;
		int i=0;
		int nseqs=0;
		for (i=0;i<ret->items;i++) {
			ret->seqsarray[i] = bseq_read(es->ks, es->opt->chunk_size, keep_comment, &ret->n_seqs[i]);
			if (es->paired == 1) ret->seqsarray2[i] = bseq_read(es->ks2, es->opt->chunk_size, keep_comment, &ret->n_seqs2[i]);
			// fprintf(stderr, "[M::%s] i= %d, read %d %d sequences\n", __func__, i, ret->n_seqs[i], ret->n_seqs2[i]);
		}
		/*fprintf(stderr,"[%s] n_seqs = %d\n",__func__,retarray[0]->n_seqs);
		fprintf(stderr,"[%s] n_seqs = %d\n",__func__,retarray[1]->n_seqs);
		fprintf(stderr,"[%s] n_seqs = %d\n",__func__,retarray[2]->n_seqs);
		fprintf(stderr,"[%s] n_seqs = %d\n",__func__,retarray[3]->n_seqs);
		fprintf(stderr,"[%s] n_seqs = %d\n",__func__,retarray[4]->n_seqs);*/
		if (ret->seqsarray[0]) return ret;
		else {
			for (i=1;i<ret->items;i++) {
				free(ret->seqsarray[i]);
				if (es->paired==1) free(ret->seqsarray2[i]);
			}
			free(ret->seqsarray);
			free(ret->n_seqs);
			if (es->paired==1) {
				free(ret->seqsarray2);
				free(ret->n_seqs2);
			}
			free(ret);
		}
	} 
	else if (step == 1) { // call kt_for
	    step_t *data = (step_t *)_data;
		kt_for(es->opt->n_threads, kt_for_cb, data, data->items);
		return data;
	} 
	else if (step == 2) { // write
		step_t *data = (step_t*)_data;
		int i,j,meepbin;
		int newlen=0;
		bseq1_t *s;
		bseq1_t *s2;
		readpair *rp;
		
		for (i = 0; i < data->items; i++) {
			s = (bseq1_t *)data->seqsarray[i];
			if(es->paired==1){
				s2 = (bseq1_t *)data->seqsarray2[i];
			}
			for (j=0; j < data->n_seqs[i];j++) {
				if (es->meeptools_mode==MMODE_SORT) {
					es->npairs+=1;
					rp=malloc(sizeof(readpair));
					rp->mee1=s[j].o_mee;
					rp->mee2=0;
					rp->rlen1=s[j].l_seq;
					rp->rlen2=0;
					rp->r1=(char *)malloc(sizeof(char) * s[j].l_seq * 4);
					rp->r2=(char *)NULL;
					copy_read(rp->r1,&s[j]);
				}
				read_stats_update(&es->in,&s[j],0);
				if (s[j].output_mode == PAIRED_OUTPUT) {
					writeReadToGzFile(es->outR1,&s[j],es->gzipped);
					read_stats_update(&es->out,&s[j],1);
				} else if (s[j].output_mode == SINGLE_OUTPUT) {
					writeReadToGzFile(es->outS,&s[j],es->gzipped);
					read_stats_update(&es->singleout,&s[j],1);
				}
				if(es->paired==1){
					if(es->meeptools_mode==MMODE_SORT) {
						rp->r2=(char *)malloc(sizeof(char) * s2[j].l_seq * 4);
						rp->mee2=s2[j].o_mee;
						rp->rlen2=s2[j].l_seq;
						copy_read(rp->r2,&s2[j]);
					}						
					read_stats_update(&es->in2,&s2[j],0);
					if (s2[j].output_mode == PAIRED_OUTPUT) {
						writeReadToGzFile(es->outR2,&s2[j],es->gzipped);
						read_stats_update(&es->out2,&s2[j],1);
					} else if (s2[j].output_mode == SINGLE_OUTPUT) {
						writeReadToGzFile(es->outS,&s2[j],es->gzipped);
						read_stats_update(&es->singleout,&s2[j],1);
					}
				}
				if (es->meeptools_mode==MMODE_SORT) {
					HASH_ADD_INT(es->pairlist,id,rp);
				}
				if ((es->meeptools_mode==MMODE_TRIM||es->meeptools_mode==MMODE_TRIM_SLOW||es->meeptools_mode==MMODE_TRIM_FAST||es->meeptools_mode==MMODE_FILTER) && (es->truncate!=0)) {
					if (es->out.nreads==es->truncate) return 0;
				}
			}
		}
		//free(data->n_seqs);free(data->seqsarray); free(data);
	}
	return 0;
}


int main(int argc, char *argv[]) {
	shared_t es;
	opt_t opt;
    char *fn;
    char *fn2;
    char *ofn;
    char *ofn2;
    char *ofns;
    int i,j,c,chunk_size;
    int meeptools_mode=0;
    int iflag=0;
    int jflag=0;
    int oflag=0;
    int pflag=0;
    int sflag=0;
    int mflag=0;
    int lflag=0;
    int tflag=0;
    int dflag=0;
    int nflag=0;
    int fflag=0;
    int zflag=0;
    int cflag=0;
    unsigned long count;
 	double mcut=100; // default mcut is 100.. fake default
 	int lcut;
 	int ncpus;
 	int truncate;
    char msg[STDSTRLEN];
    int narg=0;
    time_t end;
    double time_elapsed;
    read_stats *rs;
 
	memset(&es, 0, sizeof(shared_t));
	time(&es.begin);
	es.truncate=0;

    if (argc < 2) return usage();
    
    if (strcmp(argv[1], "append") == 0) meeptools_mode = MMODE_APPEND;
    else if (strcmp(argv[1], "filter") == 0) meeptools_mode = MMODE_FILTER;
    else if (strcmp(argv[1], "sort") == 0) meeptools_mode = MMODE_SORT;
    else if (strcmp(argv[1], "stats") == 0) meeptools_mode = MMODE_STATS;
    else if (strcmp(argv[1], "trim") == 0) meeptools_mode = MMODE_TRIM;
    else {
            fprintf(stderr, "[%s] unrecognized command '%s'\n", __func__,argv[1]);
            return 1;
    }

    if (argc == 2) return subcommand_usage(meeptools_mode);
    narg+=2; // name of the executable + subcommand

     while ((c = getopt(argc-1, argv+1, "i:j:o:p:s:m:l:t:n:fdc:zh?")) != -1)
    {
        switch (c)
        {
        case 'i':
            iflag = 1;
            fn = strdup(optarg);
            narg+=2;
            break;
        case 'j':
            jflag = 1;
            fn2 = strdup(optarg);
            narg+=2;
            break;
        case 'o':
            oflag = 1;
            ofn = strdup(optarg);
            narg+=2;
            break;
        case 'p':
            pflag = 1;
            ofn2 = strdup(optarg);
            narg+=2;
            break;
        case 's':
            sflag = 1;
            ofns = strdup(optarg);
            narg+=2;
            break;
        case 'm':
            mflag = 1;
            mcut = strtod(optarg,(char **)&msg);
            narg+=2;
            break;
 	case 'l':
            lflag = 1;
            lcut = atoi(optarg);
            narg+=2;
            break;
  	case 'n':
 	    nflag = 1;
            ncpus = atoi(optarg);
            narg+=2;
            break;
   	case 'c':
 	    cflag = 1;
            chunk_size = atoi(optarg);
            narg+=2;
            break;
	case 't':
 	    tflag = 1;
            truncate = atoi(optarg);
            narg+=2;
            break;
        case 'd':
            dflag=1;
            narg+=1;
            break;
        case 'f':
            fflag=1;
            narg+=1;
            break;
        case 'z':
            zflag=1;
            narg+=1;
            break;
        case 'h':
            narg+=1;
            return subcommand_usage(meeptools_mode);
        default:
            return subcommand_usage(meeptools_mode);
        }
    }


	if (iflag==0) return subcommand_usage(meeptools_mode);
	if ((oflag==0) && (meeptools_mode != MMODE_STATS)) return subcommand_usage(meeptools_mode);
	if ((jflag + pflag == 1) && (meeptools_mode != MMODE_STATS)) return subcommand_usage(meeptools_mode); // only jflag or pflag specified ... both needed
        if (meeptools_mode==MMODE_FILTER || meeptools_mode==MMODE_TRIM) 
        {  
    	    if (mflag==0) return subcommand_usage(meeptools_mode);
	    if (lflag==0) lcut=DEFAULTLCUT;
	}
	if (meeptools_mode==MMODE_STATS && (mflag+tflag+dflag+pflag+fflag)!=0) return subcommand_usage(meeptools_mode);
	if (meeptools_mode==MMODE_APPEND && (mflag+tflag+dflag+pflag+fflag)!=0) return subcommand_usage(meeptools_mode);
	if (tflag) es.truncate=truncate;

 	
	opt_init(&opt);
	if (nflag==1) opt.n_threads=ncpus;
	if (cflag==1) opt.chunk_size=chunk_size;
	es.opt = &opt;
	es.ks = bseq_open(fn);
	if (jflag==1) {
		es.ks2 = bseq_open(fn2);
		es.paired=1;
	} else {
		es.ks2 = (bseq_file_t *)NULL;
		es.paired=0;
	}
	read_stats_init(&es.in);
	read_stats_init(&es.in2);
	read_stats_init(&es.out);
	read_stats_init(&es.out2);
	read_stats_init(&es.singleout);
	
	if (meeptools_mode==MMODE_TRIM && dflag==1 && fflag==1) {
		subcommand_usage(MMODE_TRIM);
	}
	if (meeptools_mode==MMODE_TRIM && dflag==1) meeptools_mode=MMODE_TRIM_SLOW;
	if (meeptools_mode==MMODE_TRIM && fflag==1) meeptools_mode=MMODE_TRIM_FAST;
	es.meeptools_mode=meeptools_mode;
	es.lcut=lcut;
	es.gzipped=zflag;
	es.outR1=(gzFile *)NULL;
	es.outR2=(gzFile *)NULL;
	es.outS=(gzFile *)NULL;	
	es.npairs=0;
	es.pairlist=(readpair *)NULL;
	if (meeptools_mode==MMODE_FILTER || meeptools_mode==MMODE_TRIM || meeptools_mode==MMODE_TRIM_FAST || meeptools_mode==MMODE_TRIM_SLOW) es.mcut=mcut; // mcut has no default!
	if (meeptools_mode!=MMODE_STATS) {
		if (zflag==1) {
			es.outR1=gzopen(ofn,"wb");
			if (jflag==1) es.outR2=gzopen(ofn2,"wb");
			if (sflag==1) es.outS=gzopen(ofns,"wb");
		} else {
			es.outR1=(gzFile *)fopen(ofn,"w");
			if (jflag==1) es.outR2=(gzFile *)fopen(ofn2,"w");
			if (sflag==1) es.outS=(gzFile *)fopen(ofns,"w");			
		}
	}

	fprintf(stdout,"#Running MEEPTOOLS version %s ....\n",MEEPTOOLS_VERSION);
	
	kt_pipeline(1, kt_pipeline_cb, &es, 3); 

	fprintf(stdout, "#\n#INPUT READ1 stats...\n");
	fprintf(stdout, "File_name = %s\n",fn);
	read_stats_print(&es.in);
	if (jflag) {
		fprintf(stdout, "#\n#INPUT READ2 stats...\n");
		fprintf(stdout, "File_name = %s\n",fn2);
		read_stats_print(&es.in2);		
	}	
	if (es.meeptools_mode==MMODE_SORT) {
		fprintf(stdout,"#SORTING READS\n",__func__);
		HASH_SORT(es.pairlist,read_pair_comparator);
		readpair *rp;
		if (oflag) {
			count=0;
			for (rp=es.pairlist;rp!= (readpair *)NULL;rp=rp->hh.next) {
				// fprintf(stderr,"%s",rp->r1);
				count+=1;
				if (zflag) {
					gzwrite(es.outR1,rp->r1,strlen(rp->r1));
					if (es.paired) {
						gzwrite(es.outR2,rp->r2,strlen(rp->r2));
					}
				} else {
					fprintf((FILE *)es.outR1,rp->r1,strlen(rp->r1));
					if (es.paired) {
						fprintf((FILE *)es.outR2,rp->r2,strlen(rp->r2));
					}
				}
				if (tflag) {
					if (count==truncate) {break;}
				}
			}
		}
	} 
	else 
	{

		if (oflag) {
			fprintf(stdout, "#\n#OUTPUT READ1 stats...\n");
			fprintf(stdout, "File_name = %s\n",ofn);
			read_stats_print(&es.out);		
		}
		if (pflag) {
			fprintf(stdout, "#\n#OUTOUT READ2 stats...\n");
			fprintf(stdout, "File_name = %s\n",ofn2);
			read_stats_print(&es.out2);		
		}
		if (sflag) {
			fprintf(stdout, "#\n#OUTOUT SINGLE END READS stats...\n");
			fprintf(stdout, "File_name = %s\n",ofns);
			read_stats_print(&es.singleout);		
		}
	}



	time(&end);
	time_elapsed=difftime(end,es.begin);
	if (time_elapsed==0) time_elapsed=1;
	rs=&es.in;
	fprintf(stdout,"#MEEPTOOLS finished in %d seconds.\n",(int)time_elapsed);
	fprintf(stdout,"#MEEPTOOLS processed about %14d reads per second!!\n",(int)(rs->nreads/time_elapsed));
	fprintf(stdout,"#          which is  about %14d bases per second!!\n",(int)(rs->nbases/time_elapsed));
	fprintf(stdout,"#Thank you for using MEEPTOOLS..... VNK\n");

// close output files
	bseq_close(es.ks);
	if (es.paired==1) bseq_close(es.ks2);
	if (meeptools_mode!=MMODE_STATS) {
		if (zflag==1) {
			gzclose(es.outR1);
			if (jflag==1) gzclose(es.outR2);
			if (sflag==1) gzclose(es.outS);
		} else {
			fclose((FILE *)es.outR1);
			if (jflag==1) fclose((FILE *)es.outR2);
			if (sflag==1) fclose((FILE *)es.outS);
		}
	}	

	return 1;

}

int writeReadToGzFile(gzFile *gzf,bseq1_t *s,int gzipped)
{
	char linebuffer[3 * MAXREADLENGTH];
	char seq[MAXREADLENGTH];
	char qual[MAXREADLENGTH];
	char commentMEEP[STDSTRLEN];
	FILE *gzf1;
	if (gzf==(gzFile *)NULL) return 1;

	strnmcpy(seq,s->seq,s->trim_start,s->trim_len);
	strnmcpy(qual,s->qual,s->trim_start,s->trim_len);
	// sprintf(linebuffer,"@%s %s\n%s\n+\n%s\n",s->name,commentMEEP,seq,qual);
	if (s->comment==(char *)NULL) {
		sprintf(commentMEEP,"AQ:%.2f:MEE:%.4f:MEEP:%.4f",s->avgq,s->mee,s->meep);
	} else {
		sprintf(commentMEEP,"%s:AQ:%.2f:MEE:%.4f:MEEP:%.4f",s->comment,s->avgq,s->mee,s->meep);
	}
	sprintf(linebuffer,"@%s %s\n%s\n+\n%s\n",s->name,commentMEEP,seq,qual);
	// fprintf(stderr,linebuffer);
	if (gzipped==1) {
		gzwrite(gzf,linebuffer,strlen(linebuffer));
	} else {
		gzf1=(FILE *)gzf;
		fprintf(gzf1,linebuffer,strlen(linebuffer));
	}
	return 1;
}

int copy_read(char *r,bseq1_t *s)
{
	char seq[MAXREADLENGTH];
	char qual[MAXREADLENGTH];
	char commentMEEP[STDSTRLEN];
	strnmcpy(seq,s->seq,s->trim_start,s->trim_len);
	strnmcpy(qual,s->qual,s->trim_start,s->trim_len);
	// sprintf(linebuffer,"@%s %s\n%s\n+\n%s\n",s->name,commentMEEP,seq,qual);
	if (s->comment==(char *)NULL) {
		sprintf(commentMEEP,"AQ:%.2f:MEE:%.4f:MEEP:%.4f",s->avgq,s->mee,s->meep);
	} else {
		sprintf(commentMEEP,"%s:AQ:%.2f:MEE:%.4f:MEEP:%.4f",s->comment,s->avgq,s->mee,s->meep);
	}
	sprintf(r,"@%s %s\n%s\n+\n%s\n",s->name,commentMEEP,seq,qual);
	return 1;
}

static int usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: meeptools (Tools to calculate maximum expected error in a FASTQ read as a percentage of read length)\n");
    fprintf(stderr, "Version: %s\n\n", MEEPTOOLS_VERSION);
    fprintf(stderr, "Usage:   meeptools <command> [options]\n\n");
    fprintf(stderr, "Command: append      append MEEP score to each read\n");
    fprintf(stderr, "         filter      filter reads based on MEEP score\n");
    fprintf(stderr, "         sort        sort reads by MEEP score\n");
    fprintf(stderr, "         stats       MEEP score based stats\n");
    fprintf(stderr, "         trim        trim reads based on MEEP score\n");
    fprintf(stderr, "\n");
    return 1;
}


static int subcommand_usage(int meeptools_mode)
{
	fprintf(stderr,"meeptools_mode=%d",meeptools_mode);
	if (meeptools_mode==MMODE_APPEND)
	{
	    fprintf(stderr, "\n");
	    fprintf(stderr, "Usage:   meeptools append [options] \n\n Appends MEEP and Average Read Quality score to comment section of read id.\n\n");
	    fprintf(stderr, "Options: -i FILE    INPUT READ1 FASTQ or FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -o FILE    OUTPUT READ1 FASTQ.GZ FILE \n");
	    fprintf(stderr, "Optional:\n");
	    fprintf(stderr, "         -j FILE    INPUT READ2 FASTQ or FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -p FILE    OUTPUT READ2 FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -n INTEGER number of threads (default = 4) \n");
	    fprintf(stderr, "         -h        help\n");
	    fprintf(stderr, "\n");
        fprintf(stderr,"Examples:\n");
        fprintf(stderr,"> meeptools append -i a.fastq.gz -o a_meep.fastq.gz\n");
        fprintf(stderr,"> meeptools append -i read1.fastq -j read2.fastq -o read1_meep.fastq.gz -p read2_meep.fastq.gz\n");
    }
	else if(meeptools_mode==MMODE_FILTER)
	{	
	    fprintf(stderr, "\n");
	    fprintf(stderr, "Usage:   meeptools filter [options] \n\n Filters reads by MEEP score.\n\n");
	    fprintf(stderr, "Options: -i FILE    INPUT READ1 FASTQ or FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -o FILE    OUTPUT READ1 FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -m FLOAT   MEEP score cut off (between 0 and 20)\n");
	    fprintf(stderr, "Optional:\n");
	    fprintf(stderr, "         -j FILE    INPUT READ2 FASTQ or FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -p FILE    OUTPUT READ2 FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -n INTEGER number of threads (default = 4) \n");
	    fprintf(stderr, "         -l INTEGER read length cut off (default = %d) \n",DEFAULTLCUT);
	    fprintf(stderr, "         -s FILE    FASTQ.GZ FILE for single read (read1 or read2) meeting MEEP score cut off\n");
	    fprintf(stderr, "         -t INTEGER Truncate output number of reads\n");
	    fprintf(stderr, "         -h       help\n");
	    fprintf(stderr, "\n");
        fprintf(stderr,"Examples:\n");
        fprintf(stderr,"> meeptools filter -i a.fastq.gz -o a_meep.fastq.gz -m 1.0\n");
        fprintf(stderr,"> meeptools filter -i read1.fastq -j read2.fastq -o read1_meep.fastq.gz -p read2_meep.fastq.gz -s singletons.fastq.gz -l 90 -m 1\n");
	}
	else if (meeptools_mode==MMODE_SORT)
	{
	    fprintf(stderr, "\n");
	    fprintf(stderr, "Usage:   meeptools sort [options] \n\n Sorts reads by increasing MEEP scores (best to worse)\n\n");
	    fprintf(stderr, "Options: -i FILE    INPUT READ1 FASTQ or FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -o FILE    OUTPUT READ1 FASTQ.GZ FILE \n");
	    fprintf(stderr, "Optional:\n");
	    fprintf(stderr, "         -j FILE    INPUT READ2 FASTQ or FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -p FILE    OUTPUT READ2 FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -n INTEGER number of threads (default = 4) \n");
	    fprintf(stderr, "         -t INTEGER Truncate output number of reads\n");
	    fprintf(stderr, "         -h        help\n");
	    fprintf(stderr, "\n");
        fprintf(stderr,"Examples:\n");
        fprintf(stderr,"> meeptools sort -i a.fastq.gz -o a_meep.fastq.gz\n");
        fprintf(stderr,"> meeptools sort -i read1.fastq -j read2.fastq -o read1_meep.fastq.gz -p read2_meep.fastq.gz -t 1000\n");
	}
	else if (meeptools_mode==MMODE_STATS)
	{
	    fprintf(stderr, "\n");
	    fprintf(stderr, "Usage:   meeptools stats [options] \n\n Generates MEEP score and other stats for fastq file.\n\n");
	    fprintf(stderr, "Options: -i FILE    INPUT READ1 FASTQ or FASTQ.GZ FILE \n");
	    fprintf(stderr, "Optional:\n");
	    fprintf(stderr, "         -j FILE    INPUT READ2 FASTQ or FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -n INTEGER number of threads (default = 4) \n");
	    fprintf(stderr, "         -h        help\n");
	    fprintf(stderr, "\n");
        fprintf(stderr,"Examples:\n");
        fprintf(stderr,"> meeptools stats -i a.fastq.gz \n");
        fprintf(stderr,"> meeptools stats -i read1.fastq -j read2.fastq\n");
	}
	else if (meeptools_mode==MMODE_TRIM)
	{
	    fprintf(stderr, "\n");
	    fprintf(stderr, "Usage:   meeptools trim [options] \n\n Trims reads from 5' and 3' ends to meet MEEP cutoff.\n\n");
	    fprintf(stderr, "Options: -i FILE    INPUT READ1 FASTQ or FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -o FILE    OUTPUT READ1 FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -m FLOAT   MEEP score cut off (between 0 and 20)\n");
	    fprintf(stderr, "Optional:\n");
	    fprintf(stderr, "         -d         Less 3' aggressive trimming (slower but may result in longer trimmed reads)\n");
	    fprintf(stderr, "         -j FILE    INPUT READ2 FASTQ or FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -p FILE    OUTPUT READ2 FASTQ.GZ FILE \n");
	    fprintf(stderr, "         -n INTEGER number of threads (default = 4) \n");
	    fprintf(stderr, "         -l INTEGER read length cut off (default = %d) \n",DEFAULTLCUT);
	    fprintf(stderr, "         -s FILE    FASTQ.GZ FILE for single read (read1 or read2) meeting MEEP score cut off\n");
	    fprintf(stderr, "         -t INTEGER Truncate output number of reads\n");
	    fprintf(stderr, "         -h       help\n");
	    fprintf(stderr, "\n");
        fprintf(stderr,"Examples:\n");
        fprintf(stderr,"> meeptools trim -i a.fastq.gz -o a_meep.fastq.gz -m 1.0\n");
        fprintf(stderr,"> meeptools trim -i read1.fastq -j read2.fastq -o read1_meep.fastq.gz -p read2_meep.fastq.gz -s singletons.fastq.gz -l 90 -m 1\n");
	}
	else
	{
		return usage();
	}
  	return 1;
}


void strnmcpy(char *s1, const char *s2, size_t m,size_t n)
{
strncpy(s1,&s2[m],n);
s1[n]='\0';		
}

char * char_repeat( int n, char c ) {
  char * dest = malloc(n+1);
  memset(dest, c, n);
  dest[n] = '\0';
  return dest;
}

void read_stats_init(read_stats *rs) {
	int i;
	rs->nreads=0;
	rs->nbases=0;
	rs->qsumsum=0;
	rs->avgqsum=0.0;
	rs->meesum=0.0;
	rs->meepsum=0.0;
	rs->minrl=10000;
	rs->maxrl=0;
	rs->minaq=50.0;
	rs->maxaq=0.0;
	rs->minmeep=10.0;
	rs->maxmeep=0.0;
	for(i=0;i<7;i++) {rs->meep_bin_counts[i]=0;}
}

int read_stats_update(read_stats *rs,bseq1_t *s,int when) {
	rs->nreads+=1;
	int meepbin;
	if (when==0) {
		rs->nbases+=s->l_seq;
		rs->avgqsum+=s->o_avgq;
		rs->qsumsum+=s->o_qsum;
		rs->meesum+=s->o_mee;
		rs->meepsum+=s->o_meep;
		meepbin=(int)(s->o_meep/0.5);
		if (s->l_seq>rs->maxrl) rs->maxrl=s->l_seq;
		if (s->l_seq<rs->minrl) rs->minrl=s->l_seq;
		if (s->o_meep>rs->maxmeep) rs->maxmeep=s->o_meep;
		if (s->o_meep<rs->minmeep) rs->minmeep=s->o_meep;
		if (s->o_avgq>rs->maxaq) rs->maxaq=s->o_avgq;
		if (s->o_avgq<rs->minaq) rs->minaq=s->o_avgq;
	} else {
		rs->nbases+=s->trim_len;
		rs->avgqsum+=s->avgq;
		rs->qsumsum+=s->qsum;
		rs->meesum+=s->mee;
		rs->meepsum+=s->meep;
		meepbin=(int)(s->meep/0.5);
		if (s->trim_len>rs->maxrl) rs->maxrl=s->trim_len;
		if (s->trim_len<rs->minrl) rs->minrl=s->trim_len;
		if (s->meep>rs->maxmeep) rs->maxmeep=s->meep;
		if (s->meep<rs->minmeep) rs->minmeep=s->meep;
		if (s->avgq>rs->maxaq) rs->maxaq=s->avgq;
		if (s->avgq<rs->minaq) rs->minaq=s->avgq;		
	}
	if (meepbin>6) meepbin=6;
	rs->meep_bin_counts[meepbin]+=1;	
	return 1;
}

int read_pair_comparator(readpair *rp,readpair *rp2){
	float m1 = (rp->mee1+rp->mee2)*100.0/(rp->rlen1+rp->rlen2);
	float m2 = (rp2->mee1+rp2->mee2)*100.0/(rp2->rlen1+rp2->rlen2);
	if (m1 < m2) {
			return -1;
		} else {
			return 1;
	}
}

void read_stats_print(read_stats *rs) {
	double bin=0.0;
	double percentage=0.0;
	int i;
	char *bar;

	fprintf(stdout, "Number_of_reads = %ld\n",rs->nreads);
	fprintf(stdout, "Number_of_bases = %ld\n",rs->nbases);
	fprintf(stdout, "MEEP1_reads = %d\n",rs->meep_bin_counts[0]+rs->meep_bin_counts[1]);
	fprintf(stdout, "MEEP2_reads = %d\n",rs->meep_bin_counts[0]+rs->meep_bin_counts[1]+rs->meep_bin_counts[2]+rs->meep_bin_counts[3]);
	fprintf(stdout, "Min_RL = %d\n",rs->minrl);
	fprintf(stdout, "Max_RL = %d\n",rs->maxrl);
	fprintf(stdout, "Avg_RL = %d\n",(int)(rs->nbases*1.0/rs->nreads));
	fprintf(stdout, "Min_Avg_RQ = %.2f\n",rs->minaq);
	fprintf(stdout, "Max_Avg_RQ = %.2f\n",rs->maxaq);
	fprintf(stdout, "Avg_RQ = %.2f\n",rs->avgqsum*1.0/rs->nreads);
	fprintf(stdout, "Avg_BQ = %.2f\n",rs->qsumsum*1.0/rs->nbases);
	fprintf(stdout, "Min_MEEP = %.4f\n",rs->minmeep);
	fprintf(stdout, "Max_MEEP = %.4f\n",rs->maxmeep);
	fprintf(stdout, "Avg_MEEP = %.4f\n",rs->meepsum*1.0/rs->nreads);
	fprintf(stdout, "Overall_MEEP = %.4f\n",rs->meesum*100.0/rs->nbases);

	fprintf(stdout,"#MEEP_distribution...\n");
	for (i=0;i<6;i++) {
		bin+=0.5;
		percentage=rs->meep_bin_counts[i]*100.0/rs->nreads;
		bar=char_repeat((int)percentage,'*');
		fprintf(stdout,"# >%.2f and <=%.2f :%16ld:%5.2f%%:%-100s|\n#\n",bin-0.5,bin,rs->meep_bin_counts[i],percentage,bar);
		free(bar);
	}
	percentage=rs->meep_bin_counts[6]*100.0/rs->nreads;
	bar=char_repeat((int)percentage,'*');
	fprintf(stdout,"# >3.00            :%16ld:%5.2f%%:%-100s|\n#\n",rs->meep_bin_counts[6],percentage,bar);
	free(bar);
}
