
        M=$(ls ../../results/droso/processing/comp*.txt | wc -l)
	MAXI=$((${M}/2))
	echo "Number of comps : ${MAXI}"
        for ((i=0; i<$MAXI; i++))
        do
		echo $i
            	python3 \
                reads_to_align.py ../../results/droso/processing/comp$i.txt \
                 ../../results/droso/processing/comp$i.fq \
                0
            mkdir -p ../../results/droso/processing/STAR_alignment_$i
            STAR --genomeDir ../../results/droso/genome \
                --runMode alignReads \
                --runThreadN 8 \
                --readFilesIn ../../results/droso/processing/comp$i.fq  \
                --outFileNamePrefix ../../results/droso/processing/STAR_alignment_$i/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMattributes Standard \
                --outFilterMultimapNmax 10000 \
                --outReadsUnmapped Fastx
            bedtools intersect -wa -a ../../data/droso/TE.gtf -b ../../results/droso/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam  -split > ../../results/droso/processing/intersectionTE$i.txt
	bedtools intersect -wa -a ../../data/droso/ref.gff -b ../../results/droso/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam  -split > ../../results/droso/processing/intersectionRef$i.txt
	bedtools intersect -wb -a ../../data/droso/TE.gtf -b ../../results/droso/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam  -split > ../../results/droso/processing/seq_intersectionTE$i.txt

	
	python3 analysis_comp.py ../../results/droso/processing/comp$i.fq ../../results/droso/processing/intersectionRef$i.txt ../../results/droso/processing/intersectionTE$i.txt ../../results/droso/processing/analysis$i.txt 
	
	python3 add_ref_TE.py ../../results/droso/processing/comp$i.txt ../../results/droso/processing/seq_intersectionTE$i.txt  

	samtools index ../../results/droso/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam 
	echo "Chr :" >> ../../results/droso/processing/analysis$i.txt
	samtools view ../../results/droso/processing/STAR_alignment_$i/Aligned.sortedByCoord.out.bam | grep "RaGOO" | awk '{print($3)}' | sort -u >> ../../results/droso/processing/analysis$i.txt

        done

	grep -L "mitochondrion_genome_RaGOO" ../../results/droso/processing/analysis*.txt | sed 's/[^0-9]*//g' | sort -n  > ../../results/droso/processing/not_mito_related.txt

	grep -L "mitochondrion_genome_RaGOO" ../../results/droso/processing/analysis*.txt | xargs grep "RaGOO" | sed 's/[^0-9]*//g' | sort -n > ../../results/droso/processing/not_mito_related_but_still_matching_chr.txt

        python3 rapportAgglo.py ../../data/droso/TE.gtf ../../results/droso/processing/intersectionTE $MAXI -target > ../../results/droso/rapportHisto.txt
   	
	
