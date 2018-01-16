# PuLSE-SimulateDataset v1.3

### Testing PuLSE on large datasets

---
Testing of the PuLSE program on a range of simulated datasets is facilitated with the PuLSE-SimulateDataset program.

Taking as input a base frequency definition string, the program builds sequences based on given occurance probabilities.  For example, a standard DNA base residue can be given a 25% chance of being a A, a 25% chance of being a T, a 25% chance of being a C, and a 25% chance of being a A.  This is encoded as the string

> 0.25,0.25,0.25,0.25

Referring to the bases, A, T, C and G respectively.  To build a standard randomised DNA triplet, we repeat the above definition string 3 times and separate each repeat by a comma:

> 0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25


The randomised stretch from a trimer library would therefore be the above repeated 3 times.

>0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25

Phage libraries sometimes use "NNK" coding, N being any DNA base, whilst K is G or T.

This NNK DNA triplet can be encoded as follows:

> 0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5

In addition to specifying base frequencies, we must specify upstream and downstream markers.

A simple trimer library with CGCAGTTACC as an upstream marker and CCGATCTCTA as a downstream marker can be specified as follows:

>CGCAGTTACC,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,CCGATCTCTA

The program is invoked on the command line and takes as arguments the base frequency definition string, followed by the number of sequence "reads" that should be generated.

Whilst generating the sequences, 15 bases will be prepended and 15 appended to the end of each sequence, then, a number of bases from the beginning and end will be deleted (based on a normal distribution of mean 5 and standard deviation of 3), this is to obtain reads of differing lengths.  1% of sequences will then be corrupted by deletion of a base from somewhere in the sequence.

Running these libraries has been instrumental in bug and error finding for PuLSE.

We now present bellow the generation of 4 libraries for testing purposes.

---

### Simulated monobody library (10 randomised positions)

Based on the DNA sequence:
ACCATCTCTGGTCTGAAACCGGGTGTTGACTACACTATCACCGTATACGCAGTTACCNNKNNKNNKNNKNNKNNKNNKNNKNNKNNKCCGATCTCTATCAACTACCGTACCTCTGCAGGTGGCAGTGGGGGTAGCGG

A fastq file may be generated reflecting the above monobody section with 10 randomised amino acid positions with the following command:
> ./pulse-simulateDataset CGCAGTTACC,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,CCGATCTCTA 1000000 > simulated-pulse-10merMonobody-CGCAGTTACCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCCGATCTCTA.fastq

It may then be evaluated with PuLSE in the standard way:
> ./pulse simulated-pulse-10merMonobody-CGCAGTTACCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCCGATCTCTA.fastq CGCAGTTACCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCCGATCTCTA UAG Q

### 7mer linear library
--------------------------------------------
Based on the DNA sequence:

>AGCGCCATGGCGNNKNNKNNKNNKNNKNNKNNKGCTGCAGGTGGC

A fastq file may be generated reflecting the above linear 7-mer library with the following command:
./pulse-simulateDataset AGCGCCATGGCG,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,GCTGCAGGTGGC 1000000 > simulated-pulse-7merLinear-AGCGCCATGGCGXXXXXXXXXXXXXXXXXXXXXGCTGCAGGTGGC.fastq

It may then be evaluated with PuLSE in the standard way:

./pulse simulated-pulse-7merLinear-AGCGCCATGGCGXXXXXXXXXXXXXXXXXXXXXGCTGCAGGTGGC.fastq AGCGCCATGGCGXXXXXXXXXXXXXXXXXXXXXGCTGCAGGTGGC UAG Q


### 12-mer linear
----------------------------
Based on the DNA sequence:

> AGCGCCATGGCGNNKNNKNNKNNKNNKNNKNNKNNKNNKNNKNNKNNKGCTGCAGGTGGC

A fastq file may be generated reflecting the above linear 12-mer library with the following command:

./pulse-simulateDataset AGCGCCATGGCG,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.5,0,0.5,GCTGCAGGTGGC 1000000 > simulated-pulse-12merLinear-AGCGCCATGGCGXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGCTGCAGGTGGC.fastq

It may then be evaluated with PuLSE in the standard way:
/pulse simulated-pulse-12merLinear-AGCGCCATGGCGXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGCTGCAGGTGGC.fastq AGCGCCATGGCGXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGCTGCAGGTGGC UAG Q
