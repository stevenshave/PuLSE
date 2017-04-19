# PuLSE (Phage Library Sequence Evaluation)

### Quality control and quantification of peptide sequences explored by phage display libraries.

---

As the complexity of phage display libraries increases, the ability to screen all theoretical library members is lost. Quality control must therefore be carried out on DNA base occurance frequency derived from a subset of library reads taken from next generation sequencing.

PuLSE is designed to read in fastq next generation sequencing data describing reads from a phage display library.  Output is positional base frequencies, as well as protein residue translation position counts and normalised frequencies.  Using this output it is easy to identify skewed libraries or positions greatly enriched for certain bases.  Output is in the form of an HTML formatted report.

## Dependencies

---

- An HTML5 compatible webbrowser to view generated reports. (Tested with Google Chrome 57.0.2987.133 64-bit and Mozilla Firefox 52.0.2 64-bit).
- To decompress the sample sequencing data, a means to decompress the xz data is required.  On linux this can be achieved with unxz and on windows 7-Zip.

#### To compile from source

- A C++14 (C++1y) compiler.  Pulse has been tested using GCC 5.4.0 and Clang 3.8.0.
- To build from the Makefile, you must have GNU Make installed (tested with GNU Make 4.1)

## Precompiled binaries

***

Precompiled binaries accompany the PuLSE distribtution package for:

- 64-bit linux (Ubuntu 16.04 and compatible).
- 32 & 64-bit Windows 7 and compatible.

|Platform|32/64 bit|Sha256 sum |
|---:   |:---|:---|
|Ububtu Linux 16.04|64-bit|`c0c2081ee8b53b3c75b052d928a812b2d821422799ede5be4a113dfa8c4d05d1`|
|Windows 7|  64-bit|`a786a247e25a9ca5047e877833d0c5b0bd8a9ad5cb9a1588c37ab9e79b16e606`|
|Windows 7|  32-bit|`1d1c23bce23c150fe0d2d7c6cd58e8981cfcd57d5e5801f8744e6af2f1e74d15`|

Precompiled binaries can be found in the 'binaries subfolder'

## Compilation

---

#### Windows

Microsoft Visual Studio 15 project and solution files are included in the PuLSE distribtution.  Alternatively, a C++14 capable compiler such as ICC or MinG can be used to compile the source code using the manual Linux compilation instructions below.

#### Linux

Pulse was developed and extensively tested on Ubuntu Linux 16.04, and as such is supplied with a configure script to build a makefile compatible with the target system.  To use the build system:

>	./configure
	make
	sudo make install
	
The final make install is not necessary if you would like to run PuLSE from the build directory, or relocate the executable youself.  Make install will copy the executable to /usr/local/bin/

Alternatively, you may manually compile PuLSE via:

##### GCC

>g++ -o pulse src/PuLSE.cpp  -Wall -O3 -std=c++1y

##### clang

>clang++ -o pulse src/PuLSE.cpp  -Wall -O3 -std=c++1y

## Usage

---

If compiling from source, it is a good idea to first run PuLSE-test before running PuLSE on your sequencing data.  This will ensure that the system is working properly. PuLSE is designed to be run on output from next generation sequencing in the fastq file format and can be invoked as follows:

`./profilePhageLibrary inFile.fastq libraryDefinition[triplet residue][...]`

#### inFile.fastq

inFile.fastq is the data output from next generation sequencing of the phage library being profiled.

#### libraryDefinition

libraryDefinition is a string of characters used by PuLSE to identify flanking DNA bases surrounding the randomised library position.  With both upstream and downstream matches made, the randomised sequence between these markers is considered a library member and included in profiling statistics.  The definition takes the form of first, DNA bases encoding the upstream forward marker, a number of X characters, and finally the downstream forward marker.  The definition is reversed and complimentary bases generated in order to deal with reverse library reads in the NGS data.  The full length of the specified (and reverse, complimentary) upstream marker is always used.  However, in the case of downstream markers, only the first 3 DNA bases are used.  The example data accompanying the PuLSE distribution uses 'CGTTGCXXXXXXXXXXXXXXXTGTGCT' as the library definition, specifying that the randomised sequence of 15 DNA bases (5 amino acids) as being flanked by CGTTGC and TGTGCT.

#### [triplet residue]  (Optional parameter)

This optional parameter allows PuLSE to operate on non-standard DNA triplet to amino acid mappings.  By default, PuLSE uses the following mapping of DNA triplets to amino acid residue single letter codes:

>UUU->F, UUC->F, UUA->L, UUG->L, CUU->L, CUC->L, CUA->L, CUG->L, AUU->I, AUC->I, AUA->I, AUG->M, GUU->V, GUC->V, GUA->V, GUG->V, UCU->S, UCC->S, UCA->S, UCG->S, AGU->S, AGC->S, CCU->P, CCC->P, CCA->P, CCG->P, ACU->T, ACC->T, ACA->T, ACG->T, GCU->A, GCC->A, GCA->A, GCG->A, UAU->Y, UAC->Y, UAA->\*, UAG->\*, UGA->\*, CAU->H, CAC->H, CAA->Q, CAG->Q, GAA->E, GAG->E, AAU->N, AAC->N, AAA->K, AAG->K, GAU->D, GAC->D, UGU->C, UGC->C, UGG->W, CGU->R, CGC->R, CGA->R, CGG->R, AGA->R, AGG->R, GGU->G, GGC->G, GGA->G, GGG->G

A custom mapping may be inserted by first specifying the DNA triplet to be modified, then the single letter amino acid code as the product of the triplet. A common option for phage display systems with nonsense surpression is the alteration of the triplet UAG->\* to UAG->Q.  The change is expressed by replacing \[triplet residue\] with:

> UAG Q

Note, that multiple changes may be made to the mappings, by continuing to specify mapping on the command line.

#### Example

The PuLSE distribution is supplied with an example dataset, containing NGS data obtained from sequencing a cyclic 5-mer phage display library containing 5 randomised amino acid positions flanked by cystine residues an expressed in a system with nonsense surpression.  This data is supplied with the PuLSE distribution and is compressed with xz compression.  Before use, it must be decompressed (on linux, you may use unxz, or on windows 7-Zip).  The PuLSE library definition for this library is as follows: CGTTGCXXXXXXXXXXXXXXXTGTGCT.  Nonsense surpression is in the form of the UAG triplet remapped to produce a glutamate residue (Q).

To run PuLSE on the included dataset, supply the library definition and remap UAG to Q, we invoke PuLSE with the following command line:

> `pulse sample-pulse-5merCyclic-CGTTGCXXXXXXXXXXXXXXXTGTGCT.fastq CGTTGCXXXXXXXXXXXXXXXTGTGCT UAG Q`

PuLSE will then output a HTML report with the name `sample-pulse-5merCyclic-CGTTGCXXXXXXXXXXXXXXXTGTGCT.html`
