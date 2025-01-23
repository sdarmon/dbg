#!/bin/bash

#SBATCH --time=01:00:00 #job time limit
#SBATCH --cpus-per-task=1 #ncpu on the same node
#SBATCH --mem=1G #memory reservation

export PATH=$PATH:~/Documents/scr/

begin=`date +%s`
source environment.sh
        M=$(ls ${BASE_DIR}/comp*.txt | wc -l)
	echo "Number of comps : ${M}"
        for ((i=0; i<$M; i++))
        do
		echo $i
python3 analysis_comp.py \
	${BASE_DIR}/comp${i}.fa \
	${BASE_DIR}/aligned_ref${i}.txt \
	${BASE_DIR}/aligned_Dfam${i}.txt \
	${BASE_DIR}/cr_analysis${i}.txt \
--short
	
done

end=`date +%s`
elapsed=`expr $end - $begin`

echo "Job finished"
echo "Time taken: ${elapsed}"

