#!/bin/bash
#SBATCH --time=23:59:00 #job time limit
#SBATCH -J Mus100Analysis  #jobname
#SBATCH --output=slurm_logs/FulAnalyMus.log
#SBATCH --error=slurm_logs/FulAnalyMus.err
#SBATCH --cpus-per-task=8 #ncpu on the same node
#SBATCH --mem=16G #memory reservation

export PATH=$PATH:~/Documents/scr/

begin=`date +%s`

ENV="/beegfs/data/sdarmon/mus/environment100.sh"
bash ini_mus.sh "${ENV}"
bash inter_mus.sh "${ENV}"
bash comp_extension.sh "${ENV}"
#
#source "${ENV}"
#
#
#
#python3 ${WORK_DIR}/add_ref_TE.py \
#        ${BASE_DIR}/all_unitigs.txt \
#        ${BASE_DIR}/aligned_ref_AS_all.txt \
#        ${BASE_DIR}/aligned_Dfam_all.txt \
#        ${BASE_DIR}/aligned_Rb_all.txt \
#        ${BASE_DIR}/aligned_inter_all.txt \
#        ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.abundance \
#        ${T}

end=`date +%s`
elapsed=`expr $end - $begin`

echo "Job finished"
echo "Time taken: ${elapsed}"

