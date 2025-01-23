 
python3 \
	reads_to_align.py $1.txt \
         $1.fa \
         0
            mkdir -p STAR_alignment_$1
            STARlong --genomeDir ../../results/droso/genome \
                --runMode alignReads \
                --runThreadN 8 \
                --readFilesIn $1.fa  \
                --outFileNamePrefix STAR_alignment_$1/ \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMattributes Standard \
                --outFilterMultimapNmax 10000 \
                --outReadsUnmapped Fasta
#samtools view -h STAR_alignment_$1/Aligned.sortedByCoord.out.bam | awk '$2 == 0 || $1 ~ /^@/' | samtools view -b -o STAR_alignment_$1/output_filtered.bam
cp STAR_alignment_$1/Aligned.sortedByCoord.out.bam STAR_alignment_$1/output_filtered.bam

bedtools intersect -wa -a ../../data/droso/TE.gtf -b STAR_alignment_$1/output_filtered.bam  -split > intersectionTE$1.txt
	bedtools intersect -wa -a ../../data/droso/ref.gff -b STAR_alignment_$1/output_filtered.bam  -split > intersectionRef$1.txt
bedtools intersect -wb -a ../../data/droso/ref.gff -b STAR_alignment_$1/output_filtered.bam  -split > seq_intersectionRef$1.txt
	bedtools intersect -wb -a ../../data/droso/TE.gtf -b STAR_alignment_$1/output_filtered.bam  -split > seq_intersectionTE$1.txt

	
	python3 analysis_comp.py $1.fa intersectionRef$1.txt intersectionTE$1.txt analysis$1.txt 
	
	python3 add_ref_TE.py $1.txt seq_intersectionTE$1.txt seq_intersectionRef$1.txt ../../data/droso/graph_ovarie_1_hc_5_ovarie_2_hc_5_k41.abundance 

	samtools index STAR_alignment_$1/output_filtered.bam 
	echo "Chr :" >> analysis$1.txt
	samtools view STAR_alignment_$1/output_filtered.bam | grep "RaGOO" | awk '{print($3)}' | sort -u >> analysis$1.txt
   	
	
