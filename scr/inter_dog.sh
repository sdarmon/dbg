#!/bin/bash

        M=$(ls ../../results/dog/processing/comp*.txt | wc -l)
	N=$(ls ../../results/dog/processing/comp*_TE.txt | wc -l)
	MAXI=$((${M}-${N}))
	echo "Number of comps : ${MAXI}"
        for ((i=0; i<$MAXI; i++))
        do
		echo $i
            	python3 \
                reads_to_align.py ../../results/dog/processing/comp$i.txt \
                 ../../results/dog/processing/comp$i.fq \
                0
            mkdir -p ../../results/dog/processing/STAR_alignment_$i
            mkdir -p ../../results/dog/processing/STAR_alignment_cons_$i

            STARlong --genomeDir ../../results/dog/genome \
                --runMode alignReads \
                --runThreadN 8 \
                --readFilesIn ../../results/dog/processing/comp$i.fq  \
                --outFileNamePrefix ../../results/dog/processing/STAR_alignment_$i/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMattributes Standard \
                --outFilterMultimapNmax 10000 \
                --outReadsUnmapped Fastx

            STARlong --genomeDir ../../results/dog/te_ref \
                --runMode alignReads \
                --runThreadN 8 \
                --readFilesIn ../../results/dog/processing/comp$i.fq  \
                --outFileNamePrefix ../../results/dog/processing/STAR_alignment_cons_$i/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMattributes Standard \
                --outFilterMultimapNmax 10000 \
                --outReadsUnmapped Fastx
            
	    samtools view -h -F 256 ../../results/dog/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam | samtools view -b -o ../../results/dog/processing/STAR_alignment_$i/output_filtered.bam
	    

	samtools view -F 256 ../../results/dog/processing/STAR_alignment_cons_$i/Aligned.sortedByCoord.out.bam  > ../../results/dog/processing/seq_intersectionCons$i.txt

	    bedtools intersect -wa -a ../../data/dog/TE.gtf -b ../../results/dog/processing/STAR_alignment_$i/output_filtered.bam -split > ../../results/dog/processing/intersectionTE$i.txt
	bedtools intersect -wa -a ../../data/dog/ref.gtf -b ../../results/dog/processing/STAR_alignment_$i/output_filtered.bam -split > ../../results/dog/processing/intersectionRef$i.txt
	bedtools intersect -wb -a ../../data/dog/TE.gtf -b ../../results/dog/processing/STAR_alignment_$i/output_filtered.bam  -split > ../../results/dog/processing/seq_intersectionTE$i.txt

bedtools intersect -wb -a ../../data/dog/ref.gtf -b ../../results/dog/processing/STAR_alignment_$i/output_filtered.bam  -split > ../../results/dog/processing/seq_intersectionRef$i.txt

	python3 analysis_comp.py ../../results/dog/processing/comp$i.fq ../../results/dog/processing/intersectionRef$i.txt ../../results/dog/processing/seq_intersectionCons$i.txt ../../results/dog/processing/analysis$i.txt --custom-sep "split('\t')[8].split(';')[0][8:-1]"

python3 analysis_comp.py \
../../results/dog/processing/comp$i.fq \
../../results/dog/processing/intersectionRef$i.txt \
../../results/dog/processing/seq_intersectionCons$i.txt \
../../results/dog/processing/cr_analysis$i.txt \
--custom-sep "split('\t')[8].split(';')[0][8:-1]" \
--short
	
	python3 add_ref_TE.py ../../results/dog/processing/comp$i.txt ../../results/dog/processing/seq_intersectionTE$i.txt ../../results/dog/processing/seq_intersectionRef$i.txt ../../results/dog/processing/seq_intersectionCons$i.txt  ../../data/dog/results/graph_hc_10_dog_1_hc_10_dog_2_k41.abundance

	samtools index ../../results/dog/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam 
	
done

        python3 rapportAgglo.py ../../data/dog/TE.gtf ../../results/dog/processing/intersectionTE $MAXI -target > ../../results/dog/rapportHisto.txt
   	
	
