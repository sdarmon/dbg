#!/bin/bash

# The goal of this script is to align and analyse the neighbours of the components
# in order to find the true insertion of genes.

#Loading the environment variables and paths
source environment.sh
source ${WORK_DIR}/venv/bin/activate
### STEP 1 : Alignment of all the unitigs ###
#
###Directories
#mkdir -p ${BASE_DIR}/STAR_alignment_all
#mkdir -p ${BASE_DIR}/STAR_alignment_cons_dfam_all
#mkdir -p ${BASE_DIR}/STAR_alignment_cons_rb_all
#
###Convert the unitigs into reads for star
#cp results/outputNodes.txt ${BASE_DIR}/all_unitigs.txt
#python3 ${WORK_DIR}/reads_to_align.py ${BASE_DIR}/all_unitigs.txt ${BASE_DIR}/all_unitigs.fq 0
#
###Align the reads with starlong
#${STAR_BIN_L} \
#    --genomeDir ${BASE_DIR}/../ref \
#    --runMode alignReads \
#    --runThreadN 8 \
#    --readFilesIn ${BASE_DIR}/all_unitigs.fq \
#    --outFileNamePrefix ${BASE_DIR}/STAR_alignment_all/ \
#    --outSAMtype BAM SortedByCoordinate \
#    --outSAMunmapped Within \
#    --outSAMattributes Standard \
#    --outFilterMultimapNmax 10000 \
#    --seedPerReadNmax 10000000 \
#    --seedMultimapNmax 10000000 \
#    --outReadsUnmapped Fastx
#
#${STAR_BIN_L} \
#    --genomeDir ${BASE_DIR}/../dfam_ref \
#    --runMode alignReads \
#    --runThreadN 8 \
#    --readFilesIn ${BASE_DIR}/all_unitigs.fq \
#    --outFileNamePrefix ${BASE_DIR}/STAR_alignment_cons_dfam_all/ \
#    --outSAMtype BAM SortedByCoordinate \
#    --outSAMunmapped Within \
#    --outSAMattributes Standard \
#    --outFilterMultimapNmax 10000 \
#    --seedPerReadNmax 10000000 \
#    --seedMultimapNmax 10000000 \
#    --outReadsUnmapped Fastx
#
#${STAR_BIN_L} \
#    --genomeDir ${BASE_DIR}/../rb_ref \
#    --runMode alignReads \
#    --runThreadN 8 \
#    --readFilesIn ${BASE_DIR}/all_unitigs.fq \
#    --outFileNamePrefix ${BASE_DIR}/STAR_alignment_cons_rb_all/ \
#    --outSAMtype BAM SortedByCoordinate \
#    --outSAMunmapped Within \
#    --outSAMattributes Standard \
#    --outFilterMultimapNmax 10000 \
#    --seedPerReadNmax 10000000 \
#    --seedMultimapNmax 10000000 \
#    --outReadsUnmapped Fastx
#
#### STEP 2 : Analysis of the alignment and annotation of the nodes ###
#
#echo "Filter the reads to only keep the best mapped ones"
###Filter the reads to only keep the best mapped ones
#python3 ${WORK_DIR}/filter_bam.py \
#    -b \
#    ${BASE_DIR}/STAR_alignment_all/Aligned.sortedByCoord.out.bam \
#    ${BASE_DIR}/STAR_alignment_all/aligned_filtered.bam
#
#${BEDTOOLS_BIN} bamtobed \
#-i ${BASE_DIR}/STAR_alignment_all/aligned_filtered.bam \
#> ${BASE_DIR}/STAR_alignment_all/aligned_filtered.bed
#
#${SAMTOOLS_BIN} view -F 4 \
#    ${BASE_DIR}/STAR_alignment_all/aligned_filtered.bam \
#    | awk '{print $20}' FS="[\t:]" \
#    > ${BASE_DIR}/STAR_alignment_all/AS.txt
#
#paste <(sed 's/\t/\€/g' ${BASE_DIR}/STAR_alignment_all/aligned_filtered.bed) \
#    ${BASE_DIR}/STAR_alignment_all/AS.txt \
#    | sort > ${BASE_DIR}/STAR_alignment_all/aligned_AS.txt
##
#python3 ${WORK_DIR}/filter_bam.py \
#    ${BASE_DIR}/STAR_alignment_cons_dfam_all/Aligned.sortedByCoord.out.bam \
#    ${BASE_DIR}/aligned_Dfam_all.txt
#
#awk '{print $1"\t"$20}' FS=["\t:"] \
#${BASE_DIR}/aligned_Dfam_all.txt \
#| sort -u -t '_' -k2 -n > ${BASE_DIR}/STAR_alignment_cons_dfam_all/AS.txt
#
#python3 ${WORK_DIR}/filter_bam.py \
#    ${BASE_DIR}/STAR_alignment_cons_rb_all/Aligned.sortedByCoord.out.bam \
#    ${BASE_DIR}/aligned_Rb_all.txt
#
#awk '{print $1"\t"$20}' FS=["\t:"] \
#    ${BASE_DIR}/aligned_Rb_all.txt \
#    | sort -u -t '_' -k2 -n > ${BASE_DIR}/STAR_alignment_cons_rb_all/AS.txt
#
#echo "Intersect the reads with the reference genes"
###Intersect the reads with the reference genes
#${BEDTOOLS_BIN} intersect \
#    -wo  \
#    -a ${BASE_DIR}/STAR_alignment_all/aligned_filtered.bed \
#    -b ${GEN_GTF} \
#    -split  > ${BASE_DIR}/aligned_ref_all.txt
#
#echo "Compute the intergenics unitigs as the difference between the unitigs and the genes"
###Compute the intergenics unitigs as the difference between the unitigs and the genes
#${BEDTOOLS_BIN} intersect \
#    -v \
#    -a ${BASE_DIR}/STAR_alignment_all/aligned_filtered.bam \
#    -b ${GEN_GTF} \
#    -split \
#    | ${SAMTOOLS_BIN} view > ${BASE_DIR}/aligned_inter_all.txt
#
#echo "Adding the AS information to the nodes"
###Join the AS information to the nodes
#sed 's/\t/\€/1; s/\t/\€/1; s/\t/\€/1; s/\t/\€/1; s/\t/\€/1' \
#    ${BASE_DIR}/aligned_ref_all.txt \
#    | sort -k1 > ${BASE_DIR}/aligned_ref_sorted_all.txt
#
#join -t $'\t'\
#    ${BASE_DIR}/aligned_ref_sorted_all.txt \
#    ${BASE_DIR}/STAR_alignment_all/aligned_AS.txt \
#    | sed 's/\€/\t/g' > ${BASE_DIR}/aligned_ref_AS_all.txt
#
#
#echo "Annotate the components nodes with the TE and ref genes informations"
###Annotate the components nodes with the TE and ref genes informations
#python3 ${WORK_DIR}/add_ref_TE.py \
#        ${BASE_DIR}/all_unitigs.txt \
#        ${BASE_DIR}/aligned_ref_AS_all.txt \
#        ${BASE_DIR}/aligned_Dfam_all.txt \
#        ${BASE_DIR}/aligned_Rb_all.txt \
#        ${BASE_DIR}/aligned_inter_all.txt \
#        ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.abundance


### STEP 3 : Analysis of the neighbours of the components ###

##The number of components to compute
MAXI=$(ls ${BASE_DIR}/comp*.txt | wc -l)
echo "Number of comps : ${MAXI}"

for ((i=0; i<$MAXI; i++))
  do
  ## Dealing with the comp$i.txt file
  echo "Computing the genes of the component  ${i}..."
  mkdir -p ${BASE_DIR}/genes_of_comp$i

  ${WORK_DIR}/gene_finder.exe \
      ${BASE_DIR}/comp$i.txt \
      ${BASE_DIR}/all_unitigs_annotated.nodes \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
      ${BASE_DIR}/genes_of_comp$i
  done
