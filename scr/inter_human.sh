#!/bin/bash

        M=$(ls ../../results/human/processing/comp*.txt | wc -l)
	N=$(ls ../../results/human/processing/comp*_TE.txt | wc -l)
	MAXI=$((${M}-${N}))
	echo "Number of comps : ${MAXI}"
        for ((i=0; i<$MAXI; i++))
        do
		echo $i
            	python3 \
                reads_to_align.py ../../results/human/processing/comp$i.txt \
                 ../../results/human/processing/comp$i.fq \
                0
            mkdir -p ../../results/human/processing/STAR_alignment_$i
            STARlong --genomeDir ../../results/human/genome \
                --runMode alignReads \
                --runThreadN 8 \
                --readFilesIn ../../results/human/processing/comp$i.fq  \
                --outFileNamePrefix ../../results/human/processing/STAR_alignment_$i/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMattributes Standard \
                --outFilterMultimapNmax 10000 \
                --outReadsUnmapped Fastx
            bedtools intersect -wa -a ../../data/human/TE.gtf -b ../../results/human/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam  -split > ../../results/human/processing/intersectionTE$i.txt
	bedtools intersect -wa -a ../../data/human/ref.gtf -b ../../results/human/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam  -split > ../../results/human/processing/intersectionRef$i.txt
	bedtools intersect -wb -a ../../data/human/TE.gtf -b ../../results/human/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam  -split > ../../results/human/processing/seq_intersectionTE$i.txt

bedtools intersect -wb -a ../../data/human/ref.gtf -b ../../results/human/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam  -split > ../../results/human/processing/seq_intersectionRef$i.txt

	python3 analysis_comp.py ../../results/human/processing/comp$i.fq ../../results/human/processing/intersectionRef$i.txt ../../results/human/processing/intersectionTE$i.txt ../../results/human/processing/analysis$i.txt --custom-sep "split('\t')[8].split(';')[0][8:-1]"

python3 analysis_comp.py ../../results/human/processing/comp$i.fq ../../results/human/processing/intersectionRef$i.txt ../../results/human/processing/intersectionTE$i.txt ../../results/human/processing/cr_analysis$i.txt --custom-sep "split('\t')[8].split(';')[0][8:-1]" --short 
	
	python3 add_ref_TE.py ../../results/human/processing/comp$i.txt ../../results/human/processing/seq_intersectionTE$i.txt ../../results/human/processing/seq_intersectionRef$i.txt ../../data/human/results/graph_hc_r3_hc_r1_k41.abundance

	samtools index ../../results/human/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam 
	
done

        python3 rapportAgglo.py ../../data/human/TE.gtf ../../results/human/processing/intersectionTE $MAXI -target > ../../results/human/rapportHisto.txt
   	
	
