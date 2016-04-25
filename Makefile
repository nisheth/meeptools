CC=		gcc

CFLAGS=		-Wall -O3 -m64

CFLAGS+=	-DDEBUG -DHASH_FUNCTION=HASH_FNV

LFLAGS=		-lpthread -lm -lz

SRCS = 		meeptools.c \
		kthread.c \
		bseq.c

OBJS=		$(SRCS:.c=.o)

PROG=		meeptools

INCLUDES=	-I.

LIBS=

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

all:$(PROG)

$(PROG):$(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(PROG) $(OBJS) $(LFLAGS) $(LIBS)

.PHONY: depend clean

clean:
	rm -f *.o $(PROG) 

depend: $(SRCS)
	makedepend $(INCLUDES) $^
	
	
	

# DO NOT DELETE

meeptools_stats.o: meeptools_stats.h utils.h common.h /usr/include/malloc.h
meeptools_stats.o: /usr/include/features.h /usr/include/sys/cdefs.h
meeptools_stats.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
meeptools_stats.o: /usr/include/gnu/stubs-64.h /usr/include/math.h
meeptools_stats.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
meeptools_stats.o: /usr/include/bits/mathcalls.h /usr/include/stdlib.h
meeptools_stats.o: /usr/include/sys/types.h /usr/include/bits/types.h
meeptools_stats.o: /usr/include/bits/typesizes.h /usr/include/time.h
meeptools_stats.o: /usr/include/endian.h /usr/include/bits/endian.h
meeptools_stats.o: /usr/include/sys/select.h /usr/include/bits/select.h
meeptools_stats.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
meeptools_stats.o: /usr/include/sys/sysmacros.h
meeptools_stats.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
meeptools_stats.o: /usr/include/stdio.h /usr/include/libio.h
meeptools_stats.o: /usr/include/_G_config.h /usr/include/wchar.h
meeptools_stats.o: /usr/include/bits/wchar.h /usr/include/gconv.h
meeptools_stats.o: /usr/include/bits/stdio_lim.h
meeptools_stats.o: /usr/include/bits/sys_errlist.h /usr/include/string.h
meeptools_stats.o: /usr/include/zlib.h /usr/include/zconf.h
meeptools_stats.o: /usr/include/unistd.h /usr/include/bits/posix_opt.h
meeptools_stats.o: /usr/include/bits/confname.h /usr/include/getopt.h
meeptools_stats.o: uthash.h /usr/include/inttypes.h /usr/include/stdint.h
meeptools_stats.o: kseq.h /usr/include/ctype.h
