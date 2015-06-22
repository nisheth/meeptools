# MEEPTOOLS
---
A Collection of tools for NGS Quality Evaluation, Filtering and Trimming based on MEEP score, which considers the maximum expected error rate of the basecalls in a read rather than looking at the average PHRED Q score of the read.

```
>meeptools

Program: meeptools (Tools to calculate maximum expected error in a FASTQ read as a percentage of read length)
Version: 0.1

Usage:   meeptools <command> [options]

Command: append      append MEEP score to FASTQ file
         filter      filter FASTQ file based on MEEP score
         sort        sort reads by MEEP score
         stats       MEEP score based stats for FASTQ file
         subset      get a subset of reads based on MEEP score
         trim        trim reads based on MEEP score

```
