#MEEPTOOLS
##Basic description
MEEPTOOLS is a collection of subprograms for Quality Evaluation, Filtering and Trimming of Next Generation Sequencing reads based on MEEP score, which considers the maximum expected error rate of the basecalls in the read rather than the traditional approach of calculating the average PHRED QScore of the read. As the PHRED QScore of a basecall is logarithmically related to the probability error in the basecalling, calculating average of QScores of all bases in a read does not make mathematical sense to us. This is the motivation behind developing MEEPTOOLS.
##Downloading MEEPTOOLS
MEEPTOOLS can be downloaded from [here](https://github.com/nisheth/meeptools/archive/v0.2rc1.tar.gz).
##Installing MEEPTOOLS
After downloading the above tar.gz file, follow these steps:
```
tar xzvf meeptools-x.xx.tar.gz
cd meeptools-x.xx
make
```
If everything goes well, this is produce a "meeptools" executable in the current folder, which is ready to use!
###Dependencies
Without these libraries, the compilation of MEEPTOOLS will fail:
+ gcc standard math library (-lm) . This is generally pre-installed on all linux like systems
+ gcc pthread library (-lpthread) . This is required for taking advantage of multi-threading and making MEEPTOOLS significantly faster.
+ gcc libz library (-lz). This is required for reading and writing from gzipped files directly.

##Running MEEPTOOLS
Running the "meeptools" executable will give access to its subprograms, namely append, filter, sort, stats and trim.
```
>  ./meeptools 

Program: meeptools (Tools to calculate maximum expected error in a FASTQ read as a percentage of read length)
Version: 0.2

Usage:   meeptools <command> [options]

Command: append      append MEEP score to each read
         filter      filter reads based on MEEP score
         sort        sort reads by MEEP score
         stats       MEEP score based stats
         trim        trim reads based on MEEP score

```
Executing "meeptools" with a subprogram name as an argument will print more information about the specific subprogram. For example,
```
>  ./meeptools trim
Usage:   meeptools trim [options] 

 Trims reads from 5' and 3' ends to meet MEEP cutoff.

Options: -i FILE    INPUT READ1 FASTQ or FASTQ.GZ FILE 
         -o FILE    OUTPUT READ1 FASTQ.GZ FILE 
         -m FLOAT   MEEP score cut off (between 0 and 20)
Optional:
         -d         Less 3' aggressive trimming (slower but may result in longer trimmed reads)
         -j FILE    INPUT READ2 FASTQ or FASTQ.GZ FILE 
         -p FILE    OUTPUT READ2 FASTQ.GZ FILE 
         -n INTEGER number of threads (default = 4) 
         -l INTEGER read length cut off (default = 35) 
         -s FILE    FASTQ.GZ FILE for single read (read1 or read2) meeting MEEP score cut off
         -t INTEGER Truncate output number of reads
         -h       help

Examples:
> meeptools trim -i a.fastq.gz -o a_meep.fastq.gz -m 1.0
> meeptools trim -i read1.fastq -j read2.fastq -o read1_meep.fastq.gz -p read2_meep.fastq.gz -s singletons.fastq.gz -l 90 -m 1
```
### append
append is used to append MEE, MEEP and traditional average read quality information to the comment section of every read in the input FASTQ file.
### filter
filter is used to filter FASTQ file(s) based on the read MEEP score. All reads above a user-specified value are NOT reported in the output. This is used to generate a subset of full-length high quality reads.
### sort
sort reads all the reads in the memory, sort them in increasing order of MEEP scores and prints them out. One can truncate the output after a certain number of reads, say 1 million. This is a quick way of generating a subset of best reads you have in the input file.
### stats
stats calculates some basic statistics about the input FASTQ file, like number of reads/bases, average read length, average/min/max read quality, average/min/max MEEP scores, etc. It is a quick way of accessing the quality of the dataset and can easily be used to say, compare quality of two FASTQ files.
### trim
trim can eithere be used to
+ primarily chop off bases from the 3' end until the MEEP score criteria is satisfied (fast trimming), or
+ search the read to find the longest continuous subread which meets the MEEP score criteria (slow trimming)



