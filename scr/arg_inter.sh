#!/bin/bash
MAXI=$(ls ../../results/droso/$1/processing/comp*.txt | wc -l)
echo $MAXI
for ((i=0;i<$(($MAXI-1));i++));
do
        echo $i    
	python3 reads_to_align.py ../../results/droso/$1/processing/comp$i.txt \
                ../../results/droso/$1/processing/comp$i.fq 0
            mkdir -p ../../results/droso/$1/processing/STAR_alignment
            STAR --genomeDir ../../results/droso/genome \
                --runMode alignReads \
                --runThreadN 8 \
                --readFilesIn ../../results/droso/$1/processing/comp$i.fq  \
                --outFileNamePrefix ../../results/droso/$1/processing/STAR_alignment/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMattributes Standard \
                --outFilterMultimapNmax 10000 \
                --outReadsUnmapped Fastx
     bedtools intersect -wa -a ../../data/droso/TE.gtf  -b ../../results/droso/$1/processing/STAR_alignment/Aligned.sortedByCoord.out.bam -split  > ../../results/droso/$1/processing/intersectionTE$i.txt


     bedtools intersect -wb -a ../../data/droso/TE.gtf  -b ../../results/droso/$1/processing/STAR_alignment/Aligned.sortedByCoord.out.bam -split  > ../../results/droso/$1/processing/seq_intersectionTE$i.txt

     python3 add_ref_TE.py ../../results/droso/$1/processing/comp$i.txt ../../results/droso/$1/processing/seq_intersectionTE$i.txt

done
