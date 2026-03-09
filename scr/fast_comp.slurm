#!/bin/bash
#SBATCH --time=1-23:59:00 #job time limit
#SBATCH -J DogAna  #jobname
#SBATCH --output=slurm_logs/FastAnalyDog.log
#SBATCH --error=slurm_logs/FastAnalyDog.err
#SBATCH --cpus-per-task=8 #ncpu on the same node
#SBATCH --mem=16G #memory reservation

export PATH=$PATH:~/Documents/scr/

begin=`date +%s`

ENV="/beegfs/data/sdarmon/dog/environment_full.sh"
source "${ENV}"

#bash ${WORK_DIR}/comp_fast.sh "${ENV}"
#bash ${WORK_DIR}/te_analysis.sh "${ENV}"
#bash ${WORK_DIR}/inter_mus.sh "${ENV}"
MAXI=0 #$(ls ${BASE_DIR}/comp*.txt | wc -l)
        for file in ${BASE_DIR}/comp*.txt; do
          MAXI=$(( $MAXI + 1 ))
        done
        echo "Number of comps : ${MAXI}"
${WORK_DIR}/gene_finder_de_novo.exe \
          ${BASE_DIR}/comp \
          ${RESULTS_DIR}/graph/outputNodes.txt \
          ${RESULTS_DIR}/graph/hc_1_hc_2_k${K}_C0.05.edges \
          ${RESULTS_DIR}/graph/hc_1_hc_2_k${K}.abundance \
          ${BASE_DIR}/genes_of_comp \
          ${K} \
          ${MAXI}
end=`date +%s`
elapsed=$(( end - begin ))

echo "Job finished"
echo "Time taken: ${elapsed}"


