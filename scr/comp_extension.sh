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
pip install pandas
pip install matplotlib
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
#${BASE_DIR}/tot_cr_ana.txt
#paste <(sed 's/\t/\€/g' ${BASE_DIR}/STAR_alignment_all/aligned_filtered.bed) \
#    ${BASE_DIR}/STAR_alignment_all/AS.txt \
#    | sort > ${${BASE_DIR}/tot_cr_ana.txtBASE_DIR}/STAR_alignment_all/aligned_AS.txt

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

  echo "Consensus sequence :" >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt

  python3 seq_consensium_of_comp.py \
      ${BASE_DIR}/comp${i}_annotated.nodes \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
      ${K} >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt

  echo "Max abundance of the component :" >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt
  awk 'BEGIN {s=0} $9 > s {s=$9} END {print s}' FS="\t" ${BASE_DIR}/comp${i}_annotated.nodes >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt

  echo "Leaves abundance sum:" >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt
  awk '{s+=$9} END {print s}' FS="\t"  ${BASE_DIR}/genes_of_comp${i}/leaves.txt >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt


  done

mkdir -p ${BASE_DIR}/genes_of_comp
./gene_finder.exe \
      ${BASE_DIR}/comp \
      ${BASE_DIR}/all_unitigs_annotated.nodes \
      ${BASE_DIR}/../graph/graph_hc_1_hc_2_k41_C0.05.edges \
      ${BASE_DIR}/genes_of_comp \
      ${K} \
      ${MAXI}

for ((i=0;i<${MAXI};i++)) do python3 analysis_comp.py ${BASE_DIR}/comp${i}.fa ${BASE_DIR}/aligned_ref_AS${i}.txt ${BASE_DIR}/aligned_Dfam${i}.txt ${BASE_DIR}/cr_ana${i}.txt --short; done
python3 CR_ana_comps.py ${BASE_DIR}/cr_ana ${MAXI} > ${BASE_DIR}/tot_cr_ana.txt
paste ${BASE_DIR}/genes_of_comp/transcript_summary_comps.txt ${BASE_DIR}/tot_cr_ana.txt > ${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated.txt



#Checking consensus sequences matching TE

: > ${BASE_DIR}/genes_of_comp/consensus_seqs.fa
for ((i=0; i<$MAXI; i++))
  do
  echo ">cons_${i}" >> ${BASE_DIR}/genes_of_comp/consensus_seqs.fa
  sed '8q;d' ${BASE_DIR}/genes_of_comp$i/gene_summary.txt >> ${BASE_DIR}/genes_of_comp/consensus_seqs.fa
  done

: > ${BASE_DIR}/genes_of_comp/consensus_seqs_ab.fa
for ((i=0; i<$MAXI; i++))
  do
  echo ">cons_${i}" >> ${BASE_DIR}/genes_of_comp/consensus_seqs_ab.fa
  awk '{print $9,$2}' FS="\t" ${BASE_DIR}/comp${i}_annotated.nodes | sort -k1,1nr | head -1 | awk '{print $2}' >> ${BASE_DIR}/genes_of_comp/consensus_seqs_ab.fa
  done


mkdir -p ${BASE_DIR}/STAR_alignment_seqs_consensus
${STAR_BIN_L} \
    --genomeDir ${BASE_DIR}/../dfam_ref \
    --runMode alignReads \
    --runThreadN 8 \
    --readFilesIn ${BASE_DIR}/genes_of_comp/consensus_seqs.fa \
    --outFileNamePrefix ${BASE_DIR}/STAR_alignment_seqs_consensus/ \
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

python3 ${WORK_DIR}/filter_bam.py \
    ${BASE_DIR}/STAR_alignment_seqs_consensus/Aligned.sortedByCoord.out.bam \
    ${BASE_DIR}/aligned_seqs_cons.txt

awk '{print $1"\t"$3"\t"$10}' FS=["\t"] \
${BASE_DIR}/aligned_seqs_cons.txt \
| sort -u -t '_' -k2 -n \
| cut -c6- > ${BASE_DIR}/STAR_alignment_seqs_consensus/seqs_cons_matching_TE.txt

awk '$2!="*"{print $1}' FS=["\t"] ${BASE_DIR}/STAR_alignment_seqs_consensus/seqs_cons_matching_TE.txt > ${BASE_DIR}/id_comp_matching_TE.txt


awk 'NR==FNR {ids[$1]; next}
{
    idx = FNR - 1
    if ($5=="TE" && (idx in ids))        print $0 "\tTP";
    else if ($5=="TE" && !(idx in ids))  print $0 "\tFN";
    else if ($5!="TE" && (idx in ids))   print $0 "\tFP";
    else                                 print $0 "\tTN";
}' FS=["\t"] ${BASE_DIR}/id_comp_matching_TE.txt \
${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated.txt > ${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated_cons.txt

awk '$6=="TP"{pos+=1} $6=="FN"{faux+=1} END {print pos, faux}' FS=["\t"]  \
${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated_cons.txt > ${BASE_DIR}/genes_of_comp/seq_cons_precision.txt


## TeTools to know the count of each TE
awk '$1 ~ /^>/ {print $0, $2}' FS=['\t'_] ${DFAM_FA} | cut -c 2- > ${SPE_DIR}/rosette.fa



${TECOUNT} -rosette ${SPE_DIR}/rosette.fa  \
    -column 2 \
    -TE_fasta ${DFAM_FA} \
    -count ${SPE_DIR}/count_TE.txt \
    -bowtie2 \
    -RNA ${READS_1} \
    -RNApair ${READS_2} \
    -insert 500