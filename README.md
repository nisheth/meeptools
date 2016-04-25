# MEEPTOOLS
---
This is an older branch... To many optimization and "from-scratch" changes made in the new master branch .
A Collection of tools for NGS Quality Evaluation, Filtering and Trimming based on MEEP score, which considers the maximum expected error rate of the basecalls in a read rather than looking at the average PHRED Q score of the read.

---

## Downloading MEEPTOOLS

MEEPTOOLS can be downloaded from [here](https://github.com/nisheth/meeptools/archive/v0.1r5.tar.gz).

## Installing MEEPTOOLS

Clicking on the link above will download file meeptools-0.1r5.tar.gz. Then run,

```
>tar xzvf meeptools-0.1r5.tar.gz
```

This will create folder meeptools-0.1r5 with all the source files in it.
Now, first change to this folder and run "make" to create meeptools executable.

```
>cd meeptools-0.1r5
>make
```

## Running MEEPTOOLS

You may want to add the folder containing the meeptools executable to your PATH variable

```
>export PATH=$PATH:path_to_meeptools_folder
```

Now that meeptools is in your PATH, you may run it by simply typing "meeptools"

```
>meeptools
Program: meeptools (Tools to calculate maximum expected error in a FASTQ read as a percentage of read length)
Version: 0.1

Usage:   meeptools <command> [options]

Command: append      append MEEP score to FASTQ file
         filter      filter FASTQ file based on MEEP score
         sort        sort reads by MEEP score
         stats       MEEP score based stats for FASTQ file
         trim        trim reads based on MEEP score

```

---
## Things you need to know about MEEPTOOLS

* MEEPTOOLS takes one or more FASTQ files as input and writes output to screen and/or disk.

* FASTQ files written to disk as output are in gzipped format for efficient disk utilization.

* MEEPTOOLS expects FASTQ files to have encoded quality with the offset of 33 (Sanger), but older offset of 64 can be entered with the –f option.

* MEEPTOOLS uses the FASTQ sequence object from the [htslib](http://www.htslib.org/) to read/store/access/edit/write read sequences in the memory (and to the disk).

* More help can be obtain about the subprograms by typing "meeptools subprogram" like:

```
> meeptools trim
Usage:   meeptools trim [options]

 Trims reads from 5' and 3' ends to meet MEEP cutoff.

Options: -i FILE         FASTQ FILE(S) ... for paired-end read1 and read2 files should be comma separated
         -o FILE         FASTQ FILE(S) ... for paired-end read1 and read2 files should be comma separated
         -c FLOAT        MEEP score cut off (between 0 and 20)
Optional:
         -f INTEGER      Offset 33(default) or 64
.        -l INTEGER      read length cut off (default 35)
         -s FILE         FASTQ FILE for single read (read1 or read2) meeting MEEP score cut off
         -t INTEGER      Truncate number of reads in output files
         -q              Append read quality to read comment
         -m              Append MEE to read comment
         -h              help
```
---

## Performance testing

MEEPTOOLS performance was tested using SRA datasets: SRR502985<sup>**</sup>(set1) and SRR1290425(set2). Both these datasets are paired-end.

| Dataset | ReadLength | stats<sup>*</sup> | filter<sup>*</sup> | trim<sup>*</sup> |
|---------|------------|-------------------|--------------------|------------------|
| [set1](http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR502985)    | 101        | 3060              | 770                | 470              |
| [set2](http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1290425)    | 251        | 1050              | 207                | 183              |

<sup>*</sup> Number indicates thousands of reads processed per minute on a single core of AMD Opteron(tm) Processor 6174 machine running 64 bit Centos 6.0 operating system.

<sup>**</sup> Only first 5 million reads of SRR502985 are used for brevity.


