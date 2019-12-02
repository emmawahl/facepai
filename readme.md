# ﻿FACEPAI (& CaPReSe)
## Fast And Consistent Environmental-DNA Processing And Identification (& Conversion and Preparation of Reference Sequences)  

Scripts for executing a filtering, clustering and identification pipeline for eDNA samples  
Emma Wahlberg, Swedish Museum of Natural History, 2019, emma.wahlberg@nrm.se.

Documentation contents:
1. Software prerequisites
2. Introduction and how it works
3. Getting BOLD reference sequences
4. Preparing a FASTA-file from BOLD or other sources for analysis using CaPReSe
5. Merging and preparing reference sequences from two sources using CaPReSe
6. Making a BLAST database
7. Configuration
8. Executing FACEPAI
9. Results
10. References

### 1. Software prerequisites
The scripts FACEPAI and CaPReSe are executed in Bash in a Linux environment. The scripts have been developed and tested on Ubuntu 18.04.2 LTS.

Additional software required:  
**fastp 0.19.6** (should be installed with link in bin to allow global access) (Chen et al., 2018)  
**vsearch 2.7.1** (Rognes et al., 2016)  
**cutadapt 1.15** (Martin, 2011)  
**swarm 2.2.2** (Mahé, 2014; Mahé, 2015)  
**blastn 2.6.0** (Camacho et al., 2009)  
**makdeblastdb 2.6.0** (Camacho et al., 2009)  

The scripts are tested with version numbers above, they may work with both earlier and later versions.

### 2. Introduction and how it works
The scripts and associated software are executed in Bash in a Linux environment. Prerequisites may vary depending on operating system environment and application versions. Please see manuals for each application for appropriate modifications to the pipeline according to individual installations and requisites. Make sure that the scripts (files ending with .sh) are made executable and may write to new files (use the command `chmod u+x *.sh` from the terminal). Parts of the script are adapted from Mahé (2018). The scripts are developed and tested in a Ubuntu Linux environment.

Example data are availabe in a separate repository (https://1drv.ms/u/s!AicMLmGiK8MpiKsBcG6VGn-GV7NxnA?e=bE7NOf) (See Wahlberg, 2019 for details). Within the repository there are also a subset of data from Sigsgaard et al. (2019), where FACEPAI is compared to the method used in that paper.

The script FACEPAI first filters the forward and reversed reads for each replicate while at the same time merging them utilizing fastp with paired ended base correction option to filter out low quality reads (standard setting is keeping reads with phred >=Q15, overlap difference limit 40%, overlap length require 10%). The individual replicates (if more than one replicate is available) are then pooled. Primers and reads below minimum length are removed using cutadapt, and the sequences are then filtered again using vsearch to remove sequences with N:s. Sequences are then dereplicated before mOTU:s are clustered using swarm, and checking for chimeras is carried out using vsearch. The resulting mOTU:s are thereafter BLAST-searched against the database file using blastn.

The most time consuming part of this process is downloading the reference database. However, once that is accomplished, the process of first preparing the database using CaPReSe and thereafter using FACEPAI for pooling, merging, filtering and identifying eDNA sequences requires a minimum of hands on involvement. Multiple samples can be processed using the same reference database. Be aware that the reference sequences are saved locally, and are not automatically updated. It is recommended to have a periodical routine to update the data sources.

The resulting table can be manually curated and checked for discrepancies and ambiguities in preferable spread sheet editor. FACEPAI will return the top 10 hits for each sequence, to aid in the evaluation of each identification. If the recommended process for reference sequences file preparation described below is followed, the results will also report country of hit specimens as well as taxonomic lineage.

To install (from terminal):
Make sure you have Git installed. Change the current working directory to the location where you want the cloned directory with the FACEPAI scripts to be made.
Type:

    $ git clone https://github.com/emmawahl/facepai.git

Change to the newly created directory. Make sure the scripts are executable by typing:

    $ chmod u+x *.sh

### 3. Getting BOLD reference sequences
The script FACEPAI is constructed to format a results table using a database file retrieved from the Barcode of Life Database (BOLD). The standard FASTA-file downloaded from BOLD will not include information about location and taxonomic lineage. Therefor it is recommended to download a TSV-file (option “Combined: TSV” at BOLD website), and thereafter convert the TSV-file to a FASTA-file. The script CaPReSe can be used to convert the TSV-file to a FASTA-file. The script will at the same time automatically filter out sequences that are not assigned to a BIN URI, to assure that only validated quality sequences are kept.

    $ /PATH_TO_SCRIPT/caprese.sh -C INPUT OUTPUT

INPUT = TSV-file from BOLD
OUTPUT = name of resulting FASTA-file

### 4. Preparing a FASTA-file from BOLD or other sources for analysis using CaPReSe
The resulting file from the previous step will need some additional preparation before it is ready for analysis with FACEPAI. If you are using a different database file than a FASTA-file from BOLD and want to use the script unmodified, you need a FASTA-file with a ID followed by a pipe sign, followed by taxon name. Any additional information should also be separated with a pipe sign. Preparing a single FASTA-file for analysis following this format can be done with CaPReSe. If you want to merge the BOLD FASTA-file with GenBank data, you may skip this step.

    $ /PATH_TO_SCRIPT/caprese.sh -P NAME_OF_SOURCE INPUTFILE

NAME_OF_SOURCE = name of the source, e.g. BOLD, GenBank or any other name.
INPUTFILE = FASTA-file to be prepared.

### 5. Merging and preparing reference sequences from two sources using CaPReSe
CaPReSe can be used to merge and prepare FASTA-files for direct use with FACEPAI, keeping only unique sequences. CaPReSe will add identifiers (name of source) making it easy to distinguish the sources in the final results. This is useful for e.g. merging BOLD and GenBank data. A FASTA-file converted from a TSV-file in step 3 is an ideal source for BOLD data. A GenBank FASTA-file needs to be prepared before acceptance. Header elements must be delimited by a pipe sign (|), with GenBank ID first, followed by taxon name, and thereafter extra information (e.g linage). For creating a GenBank FASTA file with optional and customizable information it is recommended to download a GB-file and use the script GenBank_to_FASTA created by McKay & Rocap (2019). The separator in the GenBank FASTA-file should be a pipe sign, as in the BOLD FASTA-file. The headers in the options.config file needs to be changed according to what information is extracted from the GenBank file (argument `-a`).

Suggestion of commands for converting GenBank GB-file to FASTA-file, including accession number and taxon name:

    $ /genbank_to_fasta.py -i seqs.gb -s whole -a 'accessions','organism' -d pipe

The merging and preparation of two FASTA-files using CaPReSe is done in one step.

    $ /PATH_TO_SCRIPT/caprese.sh -M NAME_OF_SOURCE1 NAME_OF_SOURCE2 INPUTFILE1 INPUTFILE2

NAME_OF_SOURCE1 = name of first source, e.g. BOLD.
NAME_OF_SOURCE2 = name of second source, e.g. GenBank.
INPUTFILE1 = First FASTA-file to be merged.
INPUTFILE2 = Second FASTA-file to be merged.

### 6. Making a BLAST database
It is highly recommended to construct a BLAST database from the reference FASTA-file, this will drastically improve performance and memory use. The command for making a BLAST database is:

    $ makeblastdb -in FASTA_FILE -title "NAME_OF_DATABASE" -dbtype nucl

FASTA_FILE = the FASTA-file containing reference sequences
NAME_OF_DATABASE = the name of the database

### 7. Configuration
Configuration is carried out by editing the variables in the file “options.config”. The variables that needs to be changed accordingly are forward primer, reverse primer, minimum length and the path to the reference sequences database. Be aware that the reverse primer must be entered as reverse complemented. The table headers may be left as default if using information from BOLD, unless the reference sequences are prepared differently than described in this document, or provided from a different source with different metadata in different order.

### 8. Executing FACEPAI
FACEPAI is executed in the Bash terminal from the folder containing the FASTQ-files with reads.

    $ /PATH_TO_SCRIPT/facepai.sh SAMPLE_NAME FORWARD_IDENTIFIER.fastq REVERSE_IDENTIFIER.fastq

PATH_TO_SCRIPT = the path to where the script is stored in the file system.
SAMPLE_NAME = the name of the sample.
FORWARD_IDENTIFIER = string that identifies the file or files containing the forward reads in FASTQ format. This may be one or multiple files, as long as the identifier string is the same.
FORWARD_IDENTIFIER = same as above but for reverse reads.

Example:

    $ /home/UserName/Scripts/facepai.sh SoilSample1A _F.fastq _R.fastq

### 9. Results
The script will produce a number of files that can be used for statistics, and one tab-delimited file containing the BLAST results. The BLAST results are by default reported with the 10 top hits, along with a unique query sequence identifier, number of sequences included in the mOTU, identity in percent, e-value, query coverage in percent, source of subject (e.g. BOLD or GenBank if using concatenated files produced in CaPReSe), subject ID, BOLD BIN URI, taxon name, GenBank ID for BOLD subjects with corresponding GenBank data, country and taxonomic lineage. This may differ if another source or preparation of reference sequences are used, and if the heading settings are changed in the configuration file.

### 10. References
Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K. & Madden, T. L. 2009. BLAST+: architecture and applications. BMC bioinformatics 10: 421.

Chen, S.,Zhou, Y., Chen, . &, Gu, J. 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 17: 884–i890, doi: 10.1093/bioinformatics/bty560.

McKay, C. & Rocap, M. 2019. Convert Genbank or EMBL files to Fasta. Available at: https://rocaplab.ocean.washington.edu/tools/genbank_to_fasta/. [Accessed 18 February 2019].

Mahé, F., Rognes, T., Quince, C., de Vargas, C. & Dunthorn M. 2014. Swarm: robust and fast clustering method for amplicon-based studies. PeerJ: 2:e593, doi: 10.7717/peerj.593.

Mahé, F., Rognes, T., Quince, C., de Vargas, C. & Dunthorn M. 2015. Swarm v2: highly-scalable and high-resolution amplicon clustering. PeerJ: 3:e1420, doi: 10.7717/peerj.1420.

Mahé, F. 2018. Fred's metabarcoding pipeline. Available at: https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline. [Accessed 11 October 2018].

Martin, M. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17: 10–12, doi: 10.14806/ej.17.1.200.

Rognes, T., Flouri, T., Nichols, B., Quince, C., Mahé, F. 2016. VSEARCH: a versatile open source tool for metagenomics. PeerJ: 4:e2584. doi: 10.7717/peerj.2584.

Sigsgaard, E E, Nielsen, I B, Carl, H, Krag, M A, Knudsen, S W, Xing, Y, Holm‑Hansen, T H, Møller, P R, Thomsen, P F. Seawater environmental DNA reflects seasonality of a coastal fish community. Marine Biology 2017; doi: 10.1007/s00227-017-3147-4.

Wahlberg, E. 2019. FACEPAI – A script for fast and consistent environmental DNA processing and identification. BMC Ecology, in review.
