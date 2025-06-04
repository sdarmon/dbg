#!/bin/bash

##Parameters
SPE_NAME="fruit fly <Drosophila melanogaster>" #Work as well with ' ' instead of '_'
SPE="dro_som" #Short name of the species
K=41 #k-mer size
D_NT=10 #Distance in nucleotides to compute the weight of the nodes
T=10 #Threshold to compute the connexe components
SAMPLE_SIZE=18000000 #Number of reads to sample
REPBASE=0 #1 if the TE library is extracted from Repbase, 0 otherwise

##Path variables
DATA_DIR="/beegfs/data/sdarmon"
SPE_DIR="${DATA_DIR}/${SPE}"
GEN_GTF="${DATA_DIR}/dro/Drosophila_melanogaster.BDGP6.46.113.gtf"
GEN_FA="${DATA_DIR}/dro/Drosophila_melanogaster.BDGP6.46.dna_sm.toplevel.fa"
READS_1="${SPE_DIR}/DmG10.1Cont_Ca1_ATCACG_L005_R1.fastq"
READS_2="${SPE_DIR}/DmG10.1Cont_Ca1_ATCACG_L005_R2.fastq"
ROSETTE="${DATA_DIR}/TE_dro/fasta_rosette_droso.txt"

# The rosette fasta file is the join of the M. Fablet Rosette file (rosette_COMPLET_V4.txt)
# with the positions of the insertion of FlyBase (FBti_ID_loc.txt)
#
# Commands used :
#
# awk '$1 ~ /FBti*/ {print $0}' rosette_COMPLET_V4.txt > FBti_ID_line.txt
# sort -k1,1 FBti_ID_loc.txt > FBti_ID_loc_sorted.txt
# sort -k1,1 FBti_ID_line.txt > FBti_ID_line_sorted.txt
# join FBti_ID_loc_sorted.txt FBti_ID_line_sorted.txt > FBti_full_sorted.txt
# awk '{print ">"$1"\t>"$2"|"$6"/"$5"/"$4"|"$1}' FBti_full_sorted.txt > FBti_full_names.txt
# awk 'NR==FNR { mapping[$1]=$2; next } { if ($1 in mapping) $1 = mapping[$1]; print }' FBti_full_names.txt fasta_COMPLET_V4.txt | sed '/^--$/d' > fasta_rosette_droso.txt


##TE libraries that will be generated based on SPE_NAME
DFAM_FA="${SPE_DIR}/TE_Dfam.fa"
RB_FA="${SPE_DIR}/TE_Rb.fa"

##General directories
DOC_DIR="/beegfs/home/sdarmon/Documents"
RESULTS_DIR="${DATA_DIR}/results/${SPE}"
BASE_DIR="${DATA_DIR}/results/${SPE}/processing"
WORK_DIR="${DOC_DIR}/dbg/scr"
REPEATMASKER_DIR="${DATA_DIR}/RepeatMasker/RepeatMasker"
LIBRARY_DIR="${REPEATMASKER_DIR}/Libraries/"  #/beegfs/data/sdarmon/RepeatMasker/RepeatMasker/Libraries


##Binaries
VENV_FAMDB="${REPEATMASKER_DIR}/venv/bin/activate"
FAMDB_BIN="${REPEATMASKER_DIR}/famdb.py"
KISSPLICE_BIN="${WORK_DIR}/kissplice-binary-ubuntu-2.6.7/bin/kissplice"
STAR_BIN="${DOC_DIR}/STAR-2.7.11b/bin/Linux_x86_64/STAR"
STAR_BIN_L="${DOC_DIR}/STAR-2.7.11b/bin/Linux_x86_64/STARlong"
SAMTOOLS_BIN="${DOC_DIR}/samtools-1.9/samtools"
BEDTOOLS_BIN="${DOC_DIR}/bedtools2/bin/bedtools"
FASTP_BIN="${DOC_DIR}/fastp"

