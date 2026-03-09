#!/bin/bash
#SBATCH --time=1-23:59:00 #job time limit
#SBATCH -J DogAnalysis  #jobname
#SBATCH --output=slurm_logs/FulAnalyDog.log
#SBATCH --error=slurm_logs/FulAnalyDog.err
#SBATCH --cpus-per-task=8 #ncpu on the same node
#SBATCH --mem=16G #memory reservation

export PATH=$PATH:~/Documents/scr/

begin=`date +%s`

ENV="/beegfs/data/sdarmon/dog/environment_full.sh"
source "${ENV}"
source ${WORK_DIR}/venv/bin/activate

#bash comp_fast.sh "${ENV}"


#bash ini_mus.sh "${ENV}"
#bash inter_mus.sh "${ENV}"
#bash comp_extension.sh "${ENV}"
MAXI=0 #$(ls ${BASE_DIR}/comp*.txt | wc -l)
for file in ${BASE_DIR}/comp*.txt; do
MAXI=$(( $MAXI + 1 ))
done

for ((i=0; i<$MAXI; i++))
do
    python3 ${WORK_DIR}/add_ref_TE.py \
            ${BASE_DIR}/comp${i}.txt \
            ${BASE_DIR}/aligned_ref_AS${i}.txt \
            ${BASE_DIR}/aligned_Dfam${i}.txt \
            ${BASE_DIR}/aligned_Rb${i}.txt \
            ${BASE_DIR}/aligned_inter${i}.txt \
            ${RESULTS_DIR}/graph/hc_1_hc_2_k${K}.abundance
done


end=`date +%s`
elapsed=`expr $end - $begin`

echo "Job finished"
echo "Time taken: ${elapsed}"

