#!/bin/bash
#SBATCH --time=10:00:00 #job time limit
#SBATCH -J AnalysisGene  #jobname
#SBATCH --output=slurm_logs/analysisGene.log
#SBATCH --error=slurm_logs/analysisGene.err
#SBATCH --cpus-per-task=1 #ncpu on the same node
#SBATCH --mem=16G #memory reservation

export PATH=$PATH:~/Documents/scr/

begin=`date +%s`

#bash ini_mus.sh "/beegfs/data/sdarmon/mus/mus_environment.sh"
#bash inter_mus.sh "/beegfs/data/sdarmon/mus/mus_environment.sh"
bash comp_extension.sh "/beegfs/data/sdarmon/mus/mus_environment.sh"
end=`date +%s`
elapsed=`expr $end - $begin`

echo "Job finished"
echo "Time taken: ${elapsed}"

