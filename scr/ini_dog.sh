#!/bin/bash
#
##Compute the HC
echo "HC of the reads ..."
python3 homomorphic_compression.py ../../data/dog/dog_1_1.fq ../../data/dog/hc_dog_1.fq 5
python3 homomorphic_compression.py ../../data/dog/dog_1_2.fq ../../data/dog/hc_dog_2.fq 5

##Compute the DGB with kissplice
echo "DGB with kissplice ..."
kissplice -r ../../data/dog/hc_dog_1.fq -r ../../data/dog/hc_dog_2.fq -k 41

##Compute the weighting
echo "Weighting of the nodes..."
./graph.exe ../../data/dog/results/graph_hc_dog_1_hc_dog_2_k41.nodes ../../data/dog/results/graph_hc_dog_1_hc_dog_2_k41_C0.05.edges 10 -k 41 -o ../../data/dog/outputNodes.txt


## Compute the connexe components
echo "Agglomeration of connexe components..."
./agglo.exe ../../data/dog/outputNodes.txt ../../data/dog/results/graph_hc_dog_1_hc_dog_2_k41_C0.05.edges -c 20 -d 10  ../../results/dog -clean ../../data/dog/results/graph_hc_dog_1_hc_dog_2_k41.abundance > ../../results/dog/rapportAgglo.txt

## Compute the intersection
bash inter_dog.sh


##Fetch the neighborhood of a comp 
#bash neigh_comp.sh ../../results/dog/genome/ ../../results/dog/te_ref/ ../../data/dog/ref.gtf ../../data/dog/TE.gtf ../../data/dog/outputNodes.txt ../../data/dog/results/graph_hc_10_dog_1_hc_10_dog_2_k41_C0.05.edges ../../data/dog/results/graph_hc_10_dog_1_hc_10_dog_2_k41.abundance ../../results/dog/processing/comp2.txt ~/Documents/graph/sineCA2 41 10 10


##Build the vizitig file for visualisation
#./visitig ...
