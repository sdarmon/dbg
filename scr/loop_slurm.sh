#!/bin/bash
#

#Loading the environment variables and paths through
#the first argument or the default environment.sh
if [[ -n "$1" && "$1" == *.sh ]]; then
    source "$1"
else
    source environment.sh
fi

source ${WORK_DIR}/venv/bin/activate


### STEP 3 : Analysis of the neighbours of the components ###

##The number of components to compute
MAXI=$(ls ${BASE_DIR}/comp*.txt | wc -l)
echo "Number of comps : ${MAXI}"

for ((i=0; i<$MAXI; i++))
  do
  ## Dealing with the comp$i.txt file
  echo "Computing the genes of the component  ${i}..."
  mkdir -p ${BASE_DIR}/genes_of_comp$i
  sbatch \
	  -J genEx$i \
	  --output=slurm_logs/loop$i.log \
	  --error=slurm_logs/loop$i.err \
	  loop.slurm ${i}
 
done
