#!/bin/bash

# The goal of this script is to align and analyse the neighbours of the components
# in order to find the true insertion of genes.

#Loading the environment variables and paths through
#the first argument or the default environment.sh
if [[ -n "$1" && "$1" == *.sh ]]; then
    source "$1"
else
    source environment.sh
fi

source ${WORK_DIR}/venv/bin/activate
## STEP 1 : Alignment of all the unitigs ###
#
###Directories
#mkdir -p ${BASE_DIR}/STAR_alignment_all
#mkdir -p ${BASE_DIR}/STAR_alignment_cons_dfam_all
#mkdir -p ${BASE_DIR}/STAR_alignment_cons_rb_all
#
##Convert the unitigs into reads for star
#cp ${RESULTS_DIR}/graph/outputNodes.txt ${BASE_DIR}/all_unitigs.txt
#python3 ${WORK_DIR}/reads_to_align.py ${BASE_DIR}/all_unitigs.txt ${BASE_DIR}/all_unitigs.fa 0
#
###Align the reads with starlong
#${STAR_BIN_L} \
#    --genomeDir ${BASE_DIR}/../ref \
#    --runMode alignReads \
#    --runThreadN 8 \
#    --readFilesIn ${BASE_DIR}/all_unitigs.fa \
#    --outFileNamePrefix ${BASE_DIR}/STAR_alignment_all/ \
#    --outSAMtype BAM SortedByCoordinate \
#    --outSAMunmapped Within \
#    --outSAMattributes Standard \
#    --outFilterMultimapNmax 10000 \
#    --winReadCoverageRelativeMin 0.33 \
#    --seedMultimapNmax 10000000 \
#    --alignWindowsPerReadNmax 100000 \
#    --seedSearchStartLmax 1000 \
#    --winAnchorMultimapNmax 50000 \
#    --outReadsUnmapped Fastx

###--seedPerReadNmax 1000000 added because star asked to
${STAR_BIN_L} \
    --genomeDir ${BASE_DIR}/../dfam_ref \
    --runMode alignReads \
    --runThreadN 8 \
    --readFilesIn ${BASE_DIR}/all_unitigs.fa \
    --outFileNamePrefix ${BASE_DIR}/STAR_alignment_cons_dfam_all/ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFilterMultimapNmax 10000 \
    --winReadCoverageRelativeMin 0.5 \
    --seedMultimapNmax 10000000 \
    --alignWindowsPerReadNmax 100000 \
    --seedSearchStartLmax 1000 \
    --winAnchorMultimapNmax 50000 \
    --outReadsUnmapped Fastx


## STEP 2 : Analysis of the alignment and annotation of the nodes ###

echo "Filter the reads to only keep the best mapped ones"
##Filter the reads to only keep the best mapped ones (we don't filter the reads for this analysis)
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

python3 ${WORK_DIR}/filter_bam.py \
    ${BASE_DIR}/STAR_alignment_cons_dfam_all/Aligned.sortedByCoord.out.bam \
    ${BASE_DIR}/aligned_Dfam_all.txt

awk '{print $1"\t"$20}' FS=["\t:"] \
${BASE_DIR}/aligned_Dfam_all.txt \
| sort -u -t '_' -k2 -n > ${BASE_DIR}/STAR_alignment_cons_dfam_all/AS.txt


echo "" > ${BASE_DIR}/aligned_Rb_all.txt

echo "Intersect the reads with the reference genes"
##Intersect the reads with the reference genes
#${BEDTOOLS_BIN} intersect \
#    -wo  \
#    -bed \
#    -a ${BASE_DIR}/STAR_alignment_all/aligned_filtered.bam \
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
#sort -t $'\t' -k16,16 -n ${BASE_DIR}/aligned_ref_AS_all.txt > ${BASE_DIR}/temp.txt
#cp ${BASE_DIR}/temp.txt ${BASE_DIR}/aligned_ref_AS_all.txt

echo "Annotate the components nodes with the TE and ref genes informations"
##Annotate the components nodes with the TE and ref genes informations
python3 ${WORK_DIR}/add_ref_TE.py \
        ${BASE_DIR}/all_unitigs.txt \
        ${BASE_DIR}/aligned_ref_AS_all.txt \
        ${BASE_DIR}/aligned_Dfam_all.txt \
        ${BASE_DIR}/aligned_Rb_all.txt \
        ${BASE_DIR}/aligned_inter_all.txt \
        ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.abundance


### STEP 3 : Analysis of the neighbours of the components ###

##The number of components to compute
MAXI=$(ls ${BASE_DIR}/comp*.txt | wc -l)
echo "Number of comps : ${MAXI}"

#touch ${BASE_DIR}/expressed_genes.nodes
#touch ${BASE_DIR}/expressed_genes.edges
for ((i=0; i<$MAXI; i++))
  do
  ## Dealing with the comp$i.txt file
  echo "Computing the genes of the component  ${i}..."
  mkdir -p ${BASE_DIR}/genes_of_comp$i

  ${WORK_DIR}/gene_finder.exe \
      ${BASE_DIR}/comp${i}.txt \
      ${BASE_DIR}/all_unitigs_annotated.nodes \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
      ${BASE_DIR}/genes_of_comp$i \
      ${K}

#  python3 ${WORK_DIR}/abundance_genes.py \
#   ${BASE_DIR}/genes_of_comp${i}/expressed_genes.txt \
#   ${BASE_DIR}/all_unitigs_annotated.nodes \
#   ${K}
#
#   ###Editing the nodes and edges files
#    awk -F'\t' -v c=${i} '{print $0 "\t" c "\t" 0}' ${BASE_DIR}/genes_of_comp${i}/expressed_genes_annotated.txt >> ${BASE_DIR}/expressed_genes.nodes
#    awk -F'\t' -v c=${i} '{print $0 "\t 100 \t " 1}' ${BASE_DIR}/genes_of_comp${i}/TE.txt >> ${BASE_DIR}/expressed_genes.nodes
#    echo "${i}" >> ${BASE_DIR}/expressed_genes.nodes
#    rm ${BASE_DIR}/expressed_genes_ge.edges
#    rm ${BASE_DIR}/expressed_genes_te.edges
#    awk -F'\t' -v c=${i} -v t=${th} ' $2 > t {print c "\t" $0 }' ${BASE_DIR}/genes_of_comp${i}/expressed_genes_annotated.txt >> ${BASE_DIR}/expressed_genes_ge.edges
#    awk -F'\t' -v c=${i} '$1 != "*" {print c "\t" $0  "\t TE \t none"  }' ${BASE_DIR}/genes_of_comp${i}/TE.txt >> ${BASE_DIR}/expressed_genes_te.edges
#    join -t $'\t' ${BASE_DIR}/expressed_genes_te.edges ${BASE_DIR}/expressed_genes_ge.edges >> ${BASE_DIR}/expressed_genes.edges

    ##Sorting
    #LC_ALL=C sort -t $'\t' -k5,5 -k2,2 -k3,3gr  ${BASE_DIR}/expressed_genes.edges | awk '{print $5 "\t" $2 "\t" $3}' > ${BASE_DIR}/expressed_genes_sorted.edges

  done

#  ##Addding tab instead of @
#  sed -i "s/@/\t/g" ${BASE_DIR}/expressed_genes.edges
#  sed -i 's/%/\t/g' ${BASE_DIR}/expressed_genes.edges

