#!/bin/bash

#SBATCH --time=1:00:00 #job time limit
#SBATCH -J $1  #jobname
#SBATCH --output=/dev/stdout
#SBATCH --error=/dev/stderr
#SBATCH --cpus-per-task=1 #ncpu on the same node
#SBATCH --mem=1G #memory reservation

export PATH=$PATH:~/Documents/scr/

begin=`date +%s`
for file in /beegfs/data/sdarmon/results/mus/processing/comp*_annoted.nodes ; do ./analysis_gene_TE ${file} results/graph_hc_mus_1_hc_mus_2_k41_C0.05.edges > ${file}.report_gene_TE ; i=$((${i}+1)) ; echo -ne ${i}'\r' ; done
end=`date +%s`
elapsed=`expr $end - $begin`

echo "Job finished"
echo "Time taken: ${elapsed}"

