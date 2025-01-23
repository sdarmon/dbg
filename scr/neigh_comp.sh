#!/bin/bash
#

# Fonction pour afficher l'aide
afficher_aide() {
    echo "Usage: $0 [-h] genome_ref_dir te_ref_dir ref.gtf te.gtf te_cons.fa  graph.nodes graph.edges graph.ab comp.txt out.txt  k t d"
    echo
    echo "This function extend a comp up to its d-th neighbors while annotating all the nodes."
    echo
    echo "Options:"
    echo "  -h    Show this command"
    echo
    echo "Arguments:"
    echo "  genome_ref_dir     Directory of the reference of the genome"
    echo "  te_ref_dir         Directory of the reference of the TE"
    echo "  ref.gtf            Reference genome gtf file"
    echo "  te.gtf             Te annotation gtf file"  
    echo "  te_cons.fasta      Te consensus sequences"
    echo "  graph.nodes        The node file of the whole graph"
    echo "  graph.edges        The edge file of the whole graph"
    echo "  graph.ab           The abudance file of the whole graph"
    echo "  comp.txt           The file of the comp"
    echo "  dir_out            Output file directory"
    echo "  k                  k-mer size"
    echo "  t                  weighting threshold"
    echo "  d                  distance around the comp"

    echo
    echo "Exemple:"
    echo "  $0 ...to be added..."
    exit 0
}

# Vérification des arguments
if [ "$#" -lt 1 ]; then
    afficher_aide
fi

# Analyse des options
while getopts ":h" option; do
    case $option in
        h)
            afficher_aide
            ;;
        \?)
            echo "Option invalide : -$OPTARG" >&13
            afficher_aide
            ;;
    esac
done

# Supprimer les options traitées des arguments
shift $((OPTIND - 1))

# Vérification des arguments restants
if [ "$#" -lt 13 ]; then
    echo "Erreur : Nombre insuffisant d'arguments."
    afficher_aide
fi

# Variables
ref_ge=$1
ref_te=$2
ref=$3
te=$4
te_cons=$5
nodes=$6
edges=$7
ab=$8
in=$9
out=${10}
k=${11}
t=${12}
d=${13}
base=$(basename "$out")
source environment.sh
source venv/bin/activate
#Computing the neighborhood
./neigh.exe ${nodes} ${edges} ${in} -o ${out}_n -c 0 -d $d

#Alignment of the unitigs
python3 reads_to_align.py ${out}_n.nodes ${out}_n.unitigs 0
mkdir -p STAR_alignment_${base}
mkdir -p STAR_alignment_cons_${base}

${STAR_BIN_L} --genomeDir ${ref_ge} \
         --runMode alignReads \
         --runThreadN 8 \
         --readFilesIn ${out}_n.unitigs \
         --outFileNamePrefix STAR_alignment_${base}/ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --outFilterMultimapNmax 10000 \
         --seedMultimapNmax 10000000 \
         --outReadsUnmapped Fastx

${STAR_BIN_L} --genomeDir ${ref_te} \
         --runMode alignReads \
         --runThreadN 8 \
         --readFilesIn ${out}_n.unitigs \
         --outFileNamePrefix STAR_alignment_cons_${base}/ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --outFilterMultimapNmax 10000 \
         --outReadsUnmapped Fastx

python3 filter_bam.py STAR_alignment_${base}/Aligned.sortedByCoord.out.bam STAR_alignment_${base}/output_filtered.bam

python3 filter_bam.py STAR_alignment_cons_${base}/Aligned.sortedByCoord.out.bam  ${out}seq_intersectionCons.txt

#bedtools intersect -wb -a ${te} -b STAR_alignment_${base}/output_filtered.bam  -split >  ${out}seq_intersectionTE.txt
echo "" >   ${out}seq_intersectionConsRb.txt

${BEDTOOLS_BIN} intersect -wb -a ${ref} -b STAR_alignment_${base}/output_filtered.bam  -split > ${out}seq_intersectionRef.txt

#Adding the informations to the file
python3 add_ref_TE.py \
	      ${out}_n.nodes \
              ${out}seq_intersectionRef.txt \
	      ${out}seq_intersectionCons.txt \
	      ${out}seq_intersectionConsRb.txt \
	      ${out}seq_intersectionConsRb.txt \
      	      ${ab}

#Collecting the TEs found in the consensus
awk 'BEGIN {FS="\t"} {print $8}' ${out}_n_annoted.nodes | sed 's/; /\n/g' | sort -u  > ${out}_n_list_TE.txt
grep -Ff ${out}_n_list_TE.txt ${te_cons} > ${out}_n_list_TE_full.txt
awk 'BEGIN {RS=">"; FS="\n"; ORS=""} NR==FNR {genes[$1]; next} $1 in genes && $1 != "" {print ">"$1"\n"; for(i=2; i<=NF; i++) print toupper($i); print "\n"}' ${out}_n_list_TE_full.txt ${te_cons} > ${out}_n_TE_cons.fasta

#Generating the dbg
awk '{print $2}' ${out}_n_annoted.nodes | awk '{print ">unitig_" NR "\n" $0}' > ${out}_n.fasta
cat ${out}_n.fasta ${out}_n_TE_cons.fasta > ${out}_dbg.fa
bcalm -in ${out}_dbg.fa -out ${out} -abundance-min 1 -kmer-size ${k}




