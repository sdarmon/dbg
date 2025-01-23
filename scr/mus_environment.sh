#!/bin/bash

##Path variables
SPE="mus"
SPE_NAME="Mus_musculus"
GEN_GTF="/beegfs/data/sdarmon/${SPE}/Mus_musculus.GRCm39.112.gtf"
GEN_FA="/beegfs/data/sdarmon/${SPE}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
DFAM_FA="/beegfs/data/sdarmon/${SPE}/TE_Dfam.fa"
RB_FA="/beegfs/data/sdarmon/${SPE}/TE_Rb.fa"
READS_1="/beegfs/data/sdarmon/${SPE}/mus_10_1.fq"
READS_2="/beegfs/data/sdarmon/${SPE}/mus_10_2.fq"

##General directories
DATA_DIR="/beegfs/data/sdarmon/${SPE}"
RESULTS_DIR="/beegfs/data/sdarmon/results/${SPE}"
SPE_DIR="/beegfs/data/sdarmon/${SPE}"
BASE_DIR="/beegfs/data/sdarmon/results/${SPE}/processing"
WORK_DIR="/beegfs/home/sdarmon/Documents/dbg/scr"
LIBRARY_DIR="/beegfs/home/sdarmon/Documents/RepeatMasker/Libraries/"
REPEATMASKER_DIR="/beegfs/home/sdarmon/Documents/RepeatMasker"

##Binaries
FAMDB_BIN="/beegfs/home/sdarmon/Documents/RepeatMasker/famdb.py"
KISSPLICE_BIN="kissplice-binary-ubuntu-2.6.7/bin/kissplice"
STAR_BIN="/beegfs/home/sdarmon/Documents/STAR-2.7.11b/bin/Linux_x86_64/STAR"
STAR_BIN_L="/beegfs/home/sdarmon/Documents/STAR-2.7.11b/bin/Linux_x86_64/STARlong"
SAMTOOLS_BIN="/beegfs/home/sdarmon/Documents/samtools-1.9/samtools"
BEDTOOLS_BIN="/beegfs/home/sdarmon/Documents/bedtools2/bin/bedtools"

