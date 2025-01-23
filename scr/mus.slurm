#!/bin/bash
#SBATCH --time=02:00:00 #job time limit
#SBATCH -J musAnalysis  #jobname
#SBATCH --output=slurm_logs/musAnaly.log
#SBATCH --error=slurm_logs/musAnaly.err
#SBATCH --cpus-per-task=1 #ncpu on the same node
#SBATCH --mem=16G #memory reservation

export PATH=$PATH:~/Documents/scr/

begin=`date +%s`

#bash ini_mus.sh "/beegfs/data/sdarmon/mus/mus_environment.sh"
#bash inter_mus.sh "/beegfs/data/sdarmon/mus/mus_environment.sh"
#bash comp_extension.sh "/beegfs/data/sdarmon/mus/mus_environment.sh"
source environment.sh


### STEP 3 : Analysis of the neighbours of the components ###

##The number of components to compute
MAXI=$(ls ${BASE_DIR}/comp*.txt | wc -l)
echo "Number of comps : ${MAXI}"

touch ${BASE_DIR}/expressed_genes.nodes
touch ${BASE_DIR}/expressed_genes.edges
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

  python3 ${WORK_DIR}/abundance_genes.py \
   ${BASE_DIR}/genes_of_comp${i}/expressed_genes.txt \
   ${BASE_DIR}/all_unitigs_annotated.nodes \
   ${K}

   ###Editing the nodes and edges files
    awk -F'\t' -v c=${i} '{print $0 "\t" c "\t" 0}' ${BASE_DIR}/genes_of_comp${i}/expressed_genes_annotated.txt >> ${BASE_DIR}/expressed_genes.nodes
    awk -F'\t' -v c=${i} '{print $0 "\t 100 \t " 1}' ${BASE_DIR}/genes_of_comp${i}/TE.txt >> ${BASE_DIR}/expressed_genes.nodes
    echo "${i}" >> ${BASE_DIR}/expressed_genes.nodes
    rm ${BASE_DIR}/expressed_genes_ge.edges
    rm ${BASE_DIR}/expressed_genes_te.edges
    awk -F'\t' -v c=${i} -v t=${th} ' $2 > t {print c "\t" $0 }' ${BASE_DIR}/genes_of_comp${i}/expressed_genes_annotated.txt >> ${BASE_DIR}/expressed_genes_ge.edges
    awk -F'\t' -v c=${i} '$1 != "*" {print c "\t" $0  "\t TE \t none"  }' ${BASE_DIR}/genes_of_comp${i}/TE.txt >> ${BASE_DIR}/expressed_genes_te.edges
    join -t $'\t' ${BASE_DIR}/expressed_genes_te.edges ${BASE_DIR}/expressed_genes_ge.edges >> ${BASE_DIR}/expressed_genes.edges

    ##Sorting
    #LC_ALL=C sort -t $'\t' -k5,5 -k2,2 -k3,3gr  ${BASE_DIR}/expressed_genes.edges | awk '{print $5 "\t" $2 "\t" $3}' > ${BASE_DIR}/expressed_genes_sorted.edges

  done

  ##Addding tab instead of @
  sed -i "s/@/\t/g" ${BASE_DIR}/expressed_genes.edges
  sed -i 's/%/\t/g' ${BASE_DIR}/expressed_genes.edges


#./neigh_comp.sh \
#       	${RESULTS_DIR}/ref/ \
#	${RESULTS_DIR}/dfam_ref/ \
#	${DATA_DIR}/mus/Mus_musculus.GRCm39.112.gtf \
#	none \
#	${DATA_DIR}/mus/TE_Dfam.fa \
#	${RESULTS_DIR}/graph/graph_hc_1_hc_2_k41.nodes \
#	${RESULTS_DIR}/graph/graph_hc_1_hc_2_k41_C0.05.edges \
#	${RESULTS_DIR}/graph/graph_hc_1_hc_2_k41.abundance \
#	${BASE_DIR}/comp0.txt results/mus92 41 15 3

end=`date +%s`
elapsed=`expr $end - $begin`

echo "Job finished"
echo "Time taken: ${elapsed}"

