#!/bin/bash
#SBATCH --time=23:59:00 #job time limit
#SBATCH -J MusAnalysis  #jobname
#SBATCH --output=slurm_logs/FulAnalyMus.log
#SBATCH --error=slurm_logs/FulAnalyMus.err
#SBATCH --cpus-per-task=8 #ncpu on the same node
#SBATCH --mem=16G #memory reservation

export PATH=$PATH:~/Documents/scr/

begin=`date +%s`

ENV="/beegfs/data/sdarmon/mus/environment.sh"
source "${ENV}"
source ${WORK_DIR}/venv/bin/activate

#bash ini_mus.sh "${ENV}"
#bash inter_mus.sh "${ENV}"
#bash comp_extension.sh "${ENV}"

MAXI=$(ls ${BASE_DIR}/comp*.txt | wc -l)
#echo "Number of comps : ${MAXI}"
#
#  mkdir -p ${BASE_DIR}/genes_of_comp
#  ./gene_finder.exe \
#        ${BASE_DIR}/comp \
#        ${BASE_DIR}/all_unitigs_annotated.nodes \
#        ${BASE_DIR}/../graph/graph_hc_1_hc_2_k41_C0.05.edges \
#        ${BASE_DIR}/genes_of_comp \
#        ${K} \
#        ${MAXI}
#
#  for ((i=0;i<${MAXI};i++))
#  do
#    python3 analysis_comp.py \
#        ${BASE_DIR}/comp${i}.fa \
#        ${BASE_DIR}/aligned_ref_AS${i}.txt \
#        ${BASE_DIR}/aligned_Dfam${i}.txt \
#        ${BASE_DIR}/cr_ana${i}.txt \
#        --short; \
#  done
#
#  python3 CR_ana_comps.py ${BASE_DIR}/cr_ana ${MAXI} > ${BASE_DIR}/tot_cr_ana.txt
#  paste ${BASE_DIR}/genes_of_comp/transcript_summary_comps.txt ${BASE_DIR}/tot_cr_ana.txt > ${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated.txt
#
#
#
#mkdir -p ${BASE_DIR}/STAR_alignment_seqs_consensus
#${STAR_BIN_L} \
#    --genomeDir ${BASE_DIR}/../dfam_ref \
#    --runMode alignReads \
#    --runThreadN 8 \
#    --readFilesIn ${BASE_DIR}/genes_of_comp/consensus_seqs_poids.fa \
#    --outFileNamePrefix ${BASE_DIR}/STAR_alignment_seqs_consensus/ \
#    --outSAMtype BAM SortedByCoordinate \
#    --outSAMunmapped Within \
#    --outSAMattributes Standard \
#    --outFilterMultimapNmax 10000 \
#    --winReadCoverageRelativeMin 0.5 \
#    --seedMultimapNmax 10000000 \
#    --alignWindowsPerReadNmax 100000 \
#    --seedSearchStartLmax 1000 \
#    --winAnchorMultimapNmax 50000 \
#    --outReadsUnmapped Fastx
#
#python3 ${WORK_DIR}/filter_bam.py \
#    ${BASE_DIR}/STAR_alignment_seqs_consensus/Aligned.sortedByCoord.out.bam \
#    ${BASE_DIR}/aligned_seqs_cons_poids.txt
#
#awk '{print $1"\t"$3"\t"$10}' FS=["\t"] \
#${BASE_DIR}/aligned_seqs_cons_poids.txt \
#| sort -u -t '_' -k2 -n \
#| cut -c6- > ${BASE_DIR}/STAR_alignment_seqs_consensus/seqs_cons_matching_TE_poids.txt
#
#awk '$2!="*"{print $1}' FS=["\t"] ${BASE_DIR}/STAR_alignment_seqs_consensus/seqs_cons_matching_TE_poids.txt > ${BASE_DIR}/id_comp_matching_TE_poids.txt
#
#
#awk 'NR==FNR {ids[$1]; next}
#{
#    idx = FNR - 1
#    if ($5=="TE" && (idx in ids))        print $0 "\tTP";
#    else if ($5=="TE" && !(idx in ids))  print $0 "\tFN";
#    else if ($5!="TE" && (idx in ids))   print $0 "\tFP";
#    else                                 print $0 "\tTN";
#}' FS=["\t"] ${BASE_DIR}/id_comp_matching_TE_poids.txt \
#${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated.txt > ${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated_cons_poids.txt
#
#awk '$6=="TP"{pos+=1} $6=="FN"{faux+=1} END {p=100*pos/(pos+faux)} END {print pos, faux, p}' FS=["\t"]  ${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated_cons_poids.txt > ${BASE_DIR}/genes_of_comp/seq_cons_recall_poids.txt
#awk '$6=="TP"{pos+=1} $6=="FP"{faux+=1} END {p=100*pos/(pos+faux)} END {print pos, faux, p}' FS=["\t"]  ${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated_cons_poids.txt > ${BASE_DIR}/genes_of_comp/seq_cons_precision_poids.txt



${TECOUNT} -rosette ${SPE_DIR}/rosette.fa  \
    -column 2 \
    -TE_fasta ${DFAM_FA}.no_tab \
    -count ${SPE_DIR}/count_TE.txt \
    -bowtie2 \
    -RNA ${READS_1} \
    -RNApair ${READS_2} \
    -insert 500

end=`date +%s`
elapsed=`expr $end - $begin`



echo "Job finished"
echo "Time taken: ${elapsed}"

