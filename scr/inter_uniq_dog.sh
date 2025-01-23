awk '{print $2}' $1.nodes | awk '{print ">unitig_" NR "\n" $0}' > $1.fa

mkdir -p STAR_alignment_$1
           
STARlong --genomeDir ../../results/dog/genome \
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

bedtools intersect -wa -a ../../data/dog/TE.gtf -b STAR_alignment_$1/output_filtered.bam  -split > intersectionTE$1.txt
	bedtools intersect -wa -a ../../data/dog/ref.gtf -b STAR_alignment_$1/output_filtered.bam  -split > intersectionRef$1.txt
bedtools intersect -wb -a ../../data/dog/ref.gtf -b STAR_alignment_$1/output_filtered.bam  -split > seq_intersectionRef$1.txt
	bedtools intersect -wb -a ../../data/dog/TE.gtf -b STAR_alignment_$1/output_filtered.bam  -split > seq_intersectionTE$1.txt

	
	python3 analysis_comp.py $1.fa intersectionRef$1.txt intersectionTE$1.txt analysis$1.txt 
	
	python3 add_ref_TE.py $1.nodes seq_intersectionTE$1.txt seq_intersectionRef$1.txt ../../data/dog/results/graph_hc_10_dog_1_hc_10_dog_2_k41.abundance 

	samtools index STAR_alignment_$1/output_filtered.bam 
	  	
	
