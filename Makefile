CC=		gcc

CFLAGS=		-g -Wall -O2 -m64 

LFLAGS=		-lm -lz

SRCS = 		utils.c \
		meeptools_stats.c \
		meeptools_append.c \
		meeptools_filter.c \
		meeptools_sort.c \
		meeptools_subset.c \
		meeptools_trim.c \
		meeptools.c

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
	rm -f *.o Makefile.bak

depend: $(SRCS)
	makedepend $(INCLUDES) $^
	
	
	

