#!/bin/bash

#Loading the environment variables and paths through
#the first argument or the default environment.sh
if [[ -n "$1" && "$1" == *.sh ]]; then
    source "$1"
else
    source environment.sh
fi


##The number of components to compute
MAXI=$(ls ${BASE_DIR}/comp*.txt | wc -l)
echo "Number of comps : ${MAXI}"


##Activation the virtual environment for python
source ${WORK_DIR}/venv/bin/activate

##Compute the analysis over every component
for ((i=0; i<$MAXI; i++))
do
    echo ${i}
    mkdir -p ${BASE_DIR}/STAR_alignment_${i}
    mkdir -p ${BASE_DIR}/STAR_alignment_cons_dfam_${i}
    mkdir -p ${BASE_DIR}/STAR_alignment_cons_rb_${i}

    ##Convert the unitigs into reads for star
    python3 ${WORK_DIR}/reads_to_align.py ${BASE_DIR}/comp${i}.txt ${BASE_DIR}/comp${i}.fa 0


    ##Align the reads with starlong
    echo "Alignment of the reads"
    ${STAR_BIN_L} \
    --genomeDir ${BASE_DIR}/../ref \
    --runMode alignReads \
    --runThreadN 8 \
    --readFilesIn ${BASE_DIR}/comp${i}.fa \
    --outFileNamePrefix ${BASE_DIR}/STAR_alignment_${i}/ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFilterMultimapNmax 10000 \
    --seedMultimapNmax 10000000 \
    --outReadsUnmapped Fastx

    ${STAR_BIN_L} \
    --genomeDir ${RESULTS_DIR}/dfam_ref \
    --runMode alignReads \
    --runThreadN 8 \
    --readFilesIn ${BASE_DIR}/comp${i}.fa \
    --outFileNamePrefix ${BASE_DIR}/STAR_alignment_cons_dfam_${i}/ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFilterMultimapNmax 10000 \
    --seedMultimapNmax 10000000 \
    --outReadsUnmapped Fastx

    ${STAR_BIN_L} \
        --genomeDir ${RESULTS_DIR}/rb_ref \
        --runMode alignReads \
        --runThreadN 8 \
        --readFilesIn ${BASE_DIR}/comp${i}.fa \
        --outFileNamePrefix ${BASE_DIR}/STAR_alignment_cons_rb_${i}/ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --outFilterMultimapNmax 10000 \
        --seedMultimapNmax 10000000 \
        --outReadsUnmapped Fastx

    echo "Filter the reads to only keep the best mapped ones"
    ##Filter the reads to only keep the best mapped ones
    python3 ${WORK_DIR}/filter_bam.py \
        -b \
        ${BASE_DIR}/STAR_alignment_${i}/Aligned.sortedByCoord.out.bam \
        ${BASE_DIR}/STAR_alignment_${i}/aligned_filtered.bam

    ${BEDTOOLS_BIN} bamtobed \
    -i ${BASE_DIR}/STAR_alignment_${i}/aligned_filtered.bam \
    > ${BASE_DIR}/STAR_alignment_${i}/aligned_filtered.bed

    ${SAMTOOLS_BIN} view -F 4 \
        ${BASE_DIR}/STAR_alignment_${i}/aligned_filtered.bam \
        | awk '{print $20}' FS="[\t:]" \
        > ${BASE_DIR}/STAR_alignment_${i}/AS.txt

    paste <(sed 's/\t/\€/g' ${BASE_DIR}/STAR_alignment_${i}/aligned_filtered.bed) \
        ${BASE_DIR}/STAR_alignment_${i}/AS.txt \
        | sort > ${BASE_DIR}/STAR_alignment_${i}/aligned_AS.txt

    python3 ${WORK_DIR}/filter_bam.py \
        ${BASE_DIR}/STAR_alignment_cons_dfam_${i}/Aligned.sortedByCoord.out.bam \
        ${BASE_DIR}/aligned_Dfam${i}.txt

    awk '{print $1"\t"$20}' FS=["\t:"] \
    ${BASE_DIR}/aligned_Dfam${i}.txt \
    | sort -u -t '_' -k2 -n > ${BASE_DIR}/STAR_alignment_cons_dfam_${i}/AS.txt

    python3 ${WORK_DIR}/filter_bam.py \
        ${BASE_DIR}/STAR_alignment_cons_rb_${i}/Aligned.sortedByCoord.out.bam \
        ${BASE_DIR}/aligned_Rb${i}.txt

    awk '{print $1"\t"$20}' FS=["\t:"] \
        ${BASE_DIR}/aligned_Rb${i}.txt \
        | sort -u -t '_' -k2 -n > ${BASE_DIR}/STAR_alignment_cons_rb_${i}/AS.txt

    echo "Intersect the reads with the reference genes"
    ##Intersect the reads with the reference genes
    ${BEDTOOLS_BIN} intersect \
        -wo  \
        -a ${BASE_DIR}/STAR_alignment_${i}/aligned_filtered.bed \
        -b ${GEN_GTF} \
        -split  > ${BASE_DIR}/aligned_ref${i}.txt

    echo "Compute the intergenics unitigs as the difference between the unitigs and the genes"
    ##Compute the intergenics unitigs as the difference between the unitigs and the genes
    ${BEDTOOLS_BIN} intersect \
        -v \
        -a ${BASE_DIR}/STAR_alignment_${i}/aligned_filtered.bam \
        -b ${GEN_GTF} \
        -split \
        | ${SAMTOOLS_BIN} view > ${BASE_DIR}/aligned_inter${i}.txt

    echo "Adding the AS information to the nodes"
    ##Join the AS information to the nodes
    sed 's/\t/\€/1; s/\t/\€/1; s/\t/\€/1; s/\t/\€/1; s/\t/\€/1' \
        ${BASE_DIR}/aligned_ref${i}.txt \
        | sort -k1 > ${BASE_DIR}/aligned_ref_sorted${i}.txt

    join -t $'\t'\
        ${BASE_DIR}/aligned_ref_sorted${i}.txt \
        ${BASE_DIR}/STAR_alignment_${i}/aligned_AS.txt \
        | sed 's/\€/\t/g' > ${BASE_DIR}/aligned_ref_AS${i}.txt


    echo "Annotate the components nodes with the TE and ref genes informations"
    ##Annotate the components nodes with the TE and ref genes informations
    python3 ${WORK_DIR}/add_ref_TE.py \
            ${BASE_DIR}/comp${i}.txt \
            ${BASE_DIR}/aligned_ref_AS${i}.txt \
            ${BASE_DIR}/aligned_Dfam${i}.txt \
            ${BASE_DIR}/aligned_Rb${i}.txt \
            ${BASE_DIR}/aligned_inter${i}.txt \
            ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.abundance


    echo "Analysis of the genes and TE"
    ##Analysis of the genes and TE
    ${WORK_DIR}/analysis_gene_TE \
            ${BASE_DIR}/comp${i}_annotated.nodes \
            ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges > \
            ${BASE_DIR}/gene_TE_analysis${i}.txt
done
