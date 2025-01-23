#!/bin/bash

echo "\n\n\nGet every name of the TE actual transcribed, to be editted for every name"
grep  "$2\\$" ../../data/droso/actual_TE.fa > ../../data/droso/exp_$1_name.txt
grep "$2\\$"  ../../data/droso/FC30.table.tsv | awk {'print $4 "\t" $5'} | sort -k2 -nr > ../../data/droso/$1/sorted_count.txt

echo "\n\n\nMaking the directories"
mkdir -p ../../data/droso/$1
mkdir -p ../../results/droso/$1

echo "\n\n\nExtracting the gtf of the associated TEs"
sed 's/\$/\\$/g' ../../data/droso/exp_$1_name.txt > ../../data/droso/exp_$1_name_pattern.txt
sed -i 's/>//g'  ../../data/droso/exp_$1_name_pattern.txt
grep -f ../../data/droso/exp_$1_name_pattern.txt ../../data/droso/TE.gtf > ../../data/droso/$1/TE_$1.gtf

echo "\n\n\nAnd the reads and the unitigs"
python3 expressed_TE.py ../../data/droso/actual_TE.fa ../../data/droso/exp_$1_name.txt ../../data/droso/exp_$1.fa
kissplice -r ../../data/droso/exp_$1.fa  -c 0 -C 0 -o ../../data/droso/$1/
./graph.exe ../../data/droso/$1/graph_exp_$1_k41.nodes ../../data/droso/$1/graph_exp_$1_k41.edges 10 -k 41 -o ../../data/droso/$1/outputNodes.txt

echo "\n\n\nIntersecting it with the reads"
#bedtools intersect -wb -a ../../data/droso/$1/TE_$1.gtf -b ../../results/droso/STAR_alignment_ref_fastp/Aligned.sortedByCoord.out.bam -split > ../../results/droso/$1/original_reads.gtf

echo "\n\n\nGetting the reads id and the sequences associated"
awk '{print $16}' < ../../results/droso/$1/original_reads.gtf > ../../results/droso/$1/reads_id.txt
sed 's/^/@/g' -i ../../results/droso/$1/reads_id.txt
grep -A1 -f ../../results/droso/$1/reads_id.txt ../../data/ovarie_1_fastp_new_sec.fq > ../../results/droso/$1/almost.fa
awk 'FNR%3' < ../../results/droso/$1/almost.fa > ../../results/droso/$1/reads.fa

echo "\n\n\nAligning both the reads and the TE gotten"
kissplice -r ../../results/droso/$1/reads.fa -o ../../data/droso/$1/
./graph.exe ../../data/droso/$1/graph_reads_k41.nodes ../../data/droso/$1/graph_reads_k41_C0.05.edges 10 -k 41 -o ../../data/droso/$1/outputNodes.txt

python3 reads_to_align.py ../../data/droso/$1/outputNodes.txt ../../data/droso/$1/reads_unitigs.fa -1
kissplice -r ../../data/droso/$1/reads_unitigs.fa -r ../../data/droso/exp_$1.fa -c 0 -C 0 -o ../../data/droso/$1/
./graph.exe ../../data/droso/$1/graph_reads_unitigs_exp_$1_k41.nodes ../../data/droso/$1/graph_reads_unitigs_exp_$1_k41.edges 10 -k 41 -o ../../data/droso/$1/outputAllNodes.txt

python3 match_nodes.py ../../data/droso/$1/outputAllNodes.txt ../../data/droso/exp_$1.fa ../../data/droso/$1/outputNodes_corrected.txt

echo "\n\n\nComputing the associated subgraph"
./neigh.exe ../../data/droso/$1/outputAllNodes.txt ../../data/droso/$1/graph_reads_unitigs_exp_$1_k41.edges ../../data/droso/$1/outputNodes_corrected.txt -o ../../data/droso/$1/aplot -d 2

echo "\n\n\nCompute the origin of those nodes (TE insersion or reads)"
grep "$2\\$"  ../../data/droso/FC30.table.tsv | awk  '$5 > 0 {print $4 "\t" $5}' | sort -k2 -nr > ../../data/droso/$1/sorted_count.txt
python3 parent_TE.py ../../data/droso/$1/aplot.nodes ../../data/droso/exp_$1.fa ../../data/droso/$1/sorted_count.txt ../../data/droso/$1/aplot_enhanced.nodes ../../data/droso/$1/outputNodes.txt ../../data/droso/$1/graph_reads_k41.abundance  

sort -n ../../data/droso/$1/aplot.edges | uniq > ../../data/droso/$1/aplot_uniq.edges

