#!/bin/bash

##Parameters
SPE_NAME="Mus_musculus" #Work as well with ' ' instead of '_'
SPE="mus" #Short name of the species
K=41 #k-mer size
D_NT=10 #Distance in nucleotides to compute the weight of the nodes
T=10 #Threshold to compute the connexe components; will be updated in the script
SAMPLE_SIZE=5000000 #Number of reads to sample, default is 5000000
REPBASE=1 #1 if the TE library is extracted from Repbase, 0 otherwise
ROSETTE=""

##Path variables
DATA_DIR="/beegfs/data/sdarmon"
SPE_DIR="${DATA_DIR}/${SPE}"
GEN_GTF="${SPE_DIR}/Mus_musculus.GRCm39.112.gtf"
GEN_FA="${SPE_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
READS_1="${SPE_DIR}/ERR2680378_1.fastq" #https://www.ncbi.nlm.nih.gov/sra/?term=ERR2680378
READS_2="${SPE_DIR}/ERR2680378_2.fastq"

## Rerun skip command : "" for no skip, "1" for skip

#Skip in ini
SKIP_READ_SAMPLING=""
SKIP_TE_LIBRARIES=""
SKIP_GEN_GRAPH=""
SKIP_STAR_GENOME="1"
SKIP_STAR_ALIGN_REF="1"
SKIP_BUILD_RUST="1"


##TE libraries that will be generated based on SPE_NAME
DFAM_FA="${SPE_DIR}/TE_Dfam.fa"
RB_FA="${SPE_DIR}/TE_Rb.fa"

##General directories
DOC_DIR="/beegfs/home/sdarmon/Documents"
RESULTS_DIR="${DATA_DIR}/results/${SPE}"
BASE_DIR="${DATA_DIR}/results/${SPE}/processing"
WORK_DIR="${DOC_DIR}/dbg/scr"
REPEATMASKER_DIR="${DATA_DIR}/RepeatMasker/RepeatMasker"
LIBRARY_DIR="${REPEATMASKER_DIR}/Libraries/"


##Binaries
VENV_FAMDB="${REPEATMASKER_DIR}/venv/bin/activate"
FAMDB_BIN="${REPEATMASKER_DIR}/famdb.py"
KISSPLICE_BIN="${WORK_DIR}/kissplice-binary-ubuntu-2.6.7/bin/kissplice"
STAR_BIN="${DOC_DIR}/STAR-2.7.11b/bin/Linux_x86_64/STAR"
STAR_BIN_L="${DOC_DIR}/STAR-2.7.11b/bin/Linux_x86_64/STARlong"
SAMTOOLS_BIN="${DOC_DIR}/samtools-1.9/samtools"
BEDTOOLS_BIN="${DOC_DIR}/bedtools2/bin/bedtools"
FASTP_BIN="${DOC_DIR}/fastp"
TECOUNT="${DOC_DIR}/TEtools/TEcount.py"
