generative inference of sequence taxonomy (GIST)
================================================

version 0.8.04

rhetorica@cs.toronto.edu

http://compsysbio.org/gist

Release notes
-------------

This version represents a transitional state toward improvements in how training is handled for extremely large datasets. During the submission process to *Bioinformatics* we found that our experimental bacterial universe training set (~8600 genomes) had serious underfitting issues as a result of the very sparse sampling required to perform regression within a reasonable memory footprint. Future releases will focus explicitly on tackling the issue of batch-loading training scores, thereby removing the last memory bottleneck in the program.

Prerequisites
-------------

### Linux

Gist requires a standard x64 GNU/Linux system. It should compile under any POSIX-compliant environment if `getprocexe()` in `main.cpp` is rewritten (to-do for version 1.0!) Development was performed using GCC version 4.6. Gist will also compile on the Intel C++ Compiler version 12.1. It has not been tested with clang/llvm, BSD, AIX, WINE, or Mac OS X.

### BWA

BWA can be downloaded from its [official site](http://bio-bwa.sourceforge.net).

If you are using BWA on a large cluster, be aware that it requires write access to its own directory in some cases and may crash unexpectedly if the drive is mounted read-only. Gist only requires that "bwa" is in the PATH.

### Building databases

You can skip this part if you have already built a genome class database and an autocross weights file, and are simply installing Gist on a cluster for crunching data. In addition to the included tools, Griebel et al.'s [Flux Simulator](http://sammeth.net/confluence/display/SIM/Home) (which is Java-based) and Python 2.7 are required for generating the synthetic metagenomes used during the neural network training process. In order to run Delin and Lincomp, you will also need some taxonomic reference tables available from the NCBI FTP site, particularly:

    [ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz)
    [ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)

Extract these two archives into a working directory called 'taxonomy'. You will need to CD into this directory whenever running Delin or Lincomp.

### Other prerequisites

The other prerequisites, [alglib](http://www.alglib.net/) and [FragGeneScan](http://omics.informatics.indiana.edu/FragGeneScan/), have been embedded directly within Gist.

Installation
------------

Navigate to `gist/src/build` and type `make gist`. Copy the resulting binary (`gist`), the FGS profiles subdirectory (`fgs_profiles/`) and the configuration file (`gist.conf`) to a location in your `PATH`, e.g.

    mkdir /usr/local/gist
    cp gist /usr/local/gist
    cp gist.conf /usr/local/gist
    cp -R fgs_profiles /usr/local/gist
    export PATH=$PATH:/usr/local/gist

On shared systems it may be advisable to replace `/usr/local` with `~/bin`.

### Setting up the tools

Gist currently ships with three utilities, Delin, Lincomp, and Genepuddle2. Delin and Genepuddle2 are used during the data preparation process to annotate genomes and create synthetic metagenomes, respectively.  

1. Building Delin:

From within `gist/tools/delin`, type:

    g++ -o delin main.cpp
    cp delin /usr/local/gist

(Assuming you made the `/usr/local/gist` directory in the previous step.)

IMPORTANT: See also the "Building databases" prerequisites, above.

2. Building Lincomp:

From within `gist/tools/lincomp/build`, type:

    make lincomp
    cp lincomp /usr/local/gist

IMPORTANT: See also the "Building databases" prerequisites, above.

3. Installing Genepuddle2:

Genepuddle2 requires Python 2.7 and Griebel et al.'s [Flux Simulator](http://sammeth.net/confluence/display/SIM/Home). To use Genepuddle2, ensure these prerequisites are installed, and then go to `gist/tools/genepuddle2` and type:

    cp * /usr/local/gist

Generating simulated data
-------------------------

Detailed information on the Genepuddle pipeline for simulating metatranscriptomic reads will be available shortly. In the meantime, you are invited to use Griebel et al.'s [Flux Simulator](http://sammeth.net/confluence/display/SIM/Home), on which Genepuddle is based.

The Genepuddle pipeline simply automates running Flux Simulator for each genome and constructs fake .GTF files to enable re-use of coding sequences and real RNAseq reads in place of well-curated genomes. Note that Flux Simulator defaults to automatic polyadenylation of transcripts, which must be disabled for correct bacterial simulation.

Using Gist
----------

Type `gist -help` for general information on how to use the program. Gist only supports FASTA reads, so use a program such as [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) to convert FASTQ sequences first.

A detailed manual will be available online shortly at [our website](compsysbio.org/gist).

Examples of useful command lines:

    gist -save -classes <dir> -ts

Generates class profiles for the cDNA library files in the directory named. You can get these files from the gzip archives at `ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/` named `cds_from_genomic.fna` or `rna_from_genomic.fna`. The `-ts` switch makes Gist assume the sequences are already in frame.

    gist -record -classes <dir> -data <file.fasta> [-b 5000] [-t 8] [-disk]

Generates raw score tables for each method using the class data and saves them to disk. This is very useful in a cluster setting, as the results can be combined to expedite further processing or permit analysis by an external tool. However, the biggest advantage is that these raw scores can be re-used for both learning and validation. If this step is skipped, score generation must be done for each part.

The `-b` switch means that the program will limit itself to processing 5000 reads at a time (or any other specified number), preventing extreme memory usage at a minimal time cost. The `-disk` switch means that only one class (reference genome) will be loaded per thread at any given time, also helping to control memory usage, and finally the `-t` switch sets the number of threads to use (default = 1).

    gist -cc <name>.mcw -classes <dir> -data <train_file.fasta> -resume [-crank -1]

Generates classifier weights and saves them to `<name.mcw>` using the previously-generated score tables (see raw score generation, above.) Use the `-crank` switch to specify the taxonomic unit that should be trained for (-1 = strain, 0 = species, 1 = genus, 2 = family, etc.). A higher `-crank` value makes the classification easier to learn, but less precise; the program will always try to guess the correct strain and report its findings accordingly (default = -1).

    gist -cross <name>.mcw -classes <dir> -data <file>.fasta [-resume] [-s2q 1] [-t 8] [-o c] [-perf]

Generates actual final classifications using the learned weights. Use the `-resume` switch if you have already generated score tables (see above). Use the `-t` switch to enable multithreading (in this example, with 8 threads.) Use the `-s2q` switch to control how many hits are returned per read (at minimum) during the final pass, and the `-o` switch to control the output format (e.g. `-o c` if you want complete taxonomy information in a tab-separated value format.) The `-perf` switch will output classification sensitivity (for training and test data), and considers a read successfully classified if the correct genus or one of its parent taxa were reported. This can be changed (from reporting on genera) using the -crank switch (explained above.)

Licence
-------

Gist, Genepuddle2, Lincomp, and Delin are provided under the GNU General Public License, version 3.0. See http://www.gnu.org/copyleft/gpl.html for more information.

Gist incorporates code from [ALGLIB](http://www.alglib.net) and [FragGeneScan](http://omics.informatics.indiana.edu/FragGeneScan/), which are also provided under the GNU GPL.
