#MEEPTOOLS
##Basic description
MEEPTOOLS is a collection of subprograms for Quality Evaluation, Filtering and Trimming of Next Generation Sequencing reads based on MEEP score, which considers the maximum expected error rate of the basecalls in the read rather than the traditional approach of calculating the average PHRED QScore of the read. As the PHRED QScore of a basecall is logarithmically related to the probability error in the basecalling, calculating average of QScores of all bases in a read does not make mathematical sense to us. This is the motivation behind developing MEEPTOOLS.
##Downloading MEEPTOOLS
MEEPTOOLS can be downloaded from [here](https://github.com/nisheth/meeptools/archive/v0.1r5.tar.gz).
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
