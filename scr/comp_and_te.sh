#!/bin/bash
#Need some comments here



#Loading the environment variables and paths through
#the first argument or the default environment.sh
if [[ -n "$1" && "$1" == *.sh ]]; then
    source "$1"
    ENV=$1
else
    source environment.sh
    ENV="environment.sh"
fi

source ${VENV_FAMDB}
pip install pandas
pip install matplotlib




#Create the general directories of that species
mkdir -p $DATA_DIR
mkdir -p $RESULTS_DIR
mkdir -p $SPE_DIR
rm -r ${BASE_DIR}/*
mkdir -p $BASE_DIR
mkdir -p $RESULTS_DIR/processing
mkdir -p $RESULTS_DIR/ref
mkdir -p $RESULTS_DIR/graph
mkdir -p $RESULTS_DIR/dfam_ref
mkdir -p $RESULTS_DIR/rb_ref

  ##FastP of the reads to remove the poly(A) tails
  echo "FastP of the reads ..."

  ${FASTP_BIN} \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --trim_poly_x \
      --thread 8 \
      --poly_x_min_len 5 \
      --html ${RESULTS_DIR}/fastp_log.html \
      --in1 ${READS_1} \
      --in2 ${READS_2} \
      --out1 ${READS_1}.fastp \
      --out2 ${READS_2}.fastp

  ##Compute the HC
  echo "HC of the reads ..."
  python3 homomorphic_compression.py  ${READS_1}.fastp ${DATA_DIR}/${SPE}/hc_1.fq 5
  python3 homomorphic_compression.py  ${READS_2}.fastp ${DATA_DIR}/${SPE}/hc_2.fq 5

if [[ ${SKIP_GEN_GRAPH} == "" ]]; then

  ##Compute the DGB with kissplice
  ulimit -s unlimited
  echo "DGB with kissplice ..."
  $KISSPLICE_BIN \
      -r ${DATA_DIR}/${SPE}/hc_1.fq \
      -r ${DATA_DIR}/${SPE}/hc_2.fq \
      -k ${K} \
      -o ${RESULTS_DIR}/graph

  ##Compute the weighting
  echo "Weighting of the nodes..."
  g++ -g graph.cpp ponderation.cpp -o graph.exe
  ${WORK_DIR}/graph.exe \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.nodes \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
      ${D_NT} \
      -k ${K}  \
      -o ${RESULTS_DIR}/graph/outputNodes.txt


  ##Compute the threshold
  echo "Threshold of the nodes..."
  T=$(python3 ${WORK_DIR}/plot.py ${RESULTS_DIR}/graph/outputNodes.txt top1)
  # Update the T variable in the environment.sh file
  sed -i "s/^T=.*/T=${T}/" "${ENV}"
  echo "T=${T}"

  ## Compute the connexe components
  echo "Agglomeration of connexe components..."
  g++ -g ${WORK_DIR}/graph.cpp ${WORK_DIR}/agglo.cpp -o ${WORK_DIR}/agglo.exe

  ${WORK_DIR}/agglo.exe \
      ${RESULTS_DIR}/graph/outputNodes.txt \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
      -c ${T} \
      -d ${D_NT} \
      ${RESULTS_DIR} \
      -clean ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.abundance > ${RESULTS_DIR}/rapportAgglo.txt
  #  agglo.exe ${RESULTS_DIR}/graph/outputNodes.txt  ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges -c ${T} -d ${D_NT} ${RESULTS_DIR} -clean ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.abundance > ${RESULTS_DIR}/rapportAgglo.txt
fi


if [[ ${SKIP_TE_LIBRARIES} == "" ]]; then
  # Extraction des TEs de la base de données Dfam
  echo "Extraction of the TEs from Dfam ..."
  ${FAMDB_BIN} -i ${LIBRARY_DIR}/famdb/ families \
  --include-class-in-name \
  --curated \
  --ancestors \
  --descendants \
  "${SPE_NAME}" --format fasta_name \
   | sed 's/#/\t/g' \
   | sed 's/ @/\t/g' \
   | sed 's/^>/>dfam_/g' > ${DFAM_FA}
fi

if [[ ${SKIP_STAR_GENOME} == "" ]]; then
  ## Compute the genome STAR ref

  $STAR_BIN \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir ${RESULTS_DIR}/dfam_ref \
  --genomeFastaFiles $DFAM_FA \
  --genomeSAindexNbases 9 \
  --genomeSAsparseD 4 \
  --genomeChrBinNbits 16
fi

if [[ ${SKIP_BUILD_RUST} == "" ]]; then
  ##Build the bin for the graph
  echo "Building the bin for the graph..."
  cargo build --release --manifest-path ${WORK_DIR}/graph/Cargo.toml
  cp ${WORK_DIR}/graph/target/release/graph ${WORK_DIR}/graph.exe

  ##Build the bin for the gene_finder function
  echo "Building the bin for the gene_finder function..."
  cargo build --release --manifest-path ${WORK_DIR}/gene_finder/Cargo.toml
  cp ${WORK_DIR}/gene_finder/target/release/gene_finder ${WORK_DIR}/gene_finder.exe
fi


chmod +x ${TECOUNT}


##The number of components to compute
MAXI=$(ls ${BASE_DIR}/comp*.txt | wc -l)
echo "Number of comps : ${MAXI}"

##Activation the virtual environment for python
source ${WORK_DIR}/venv/bin/activate

##Compute the analysis over every component
for ((i=0; i<$MAXI; i++))
do
    echo ${i}

    ${STAR_BIN_L} \
    --genomeDir ${RESULTS_DIR}/dfam_ref \
    --runMode alignReads \
    --runThreadN 8 \
    --readFilesIn ${BASE_DIR}/comp${i}.fa \
    --outFileNamePrefix ${BASE_DIR}/STAR_alignment_cons_dfam_${i}/ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --seedMultimapNmax 10000000 \
    --outReadsUnmapped Fastx \
    --outFilterMultimapNmax 10000 \
    --winReadCoverageRelativeMin 0.5 \
    --alignWindowsPerReadNmax 50000 \
    --seedSearchStartLmax 1000 \
    --winAnchorMultimapNmax 50000

    python3 ${WORK_DIR}/filter_bam.py \
        ${BASE_DIR}/STAR_alignment_cons_dfam_${i}/Aligned.sortedByCoord.out.bam \
        ${BASE_DIR}/aligned_Dfam${i}.txt

    awk '{print $1"\t"$20}' FS=["\t:"] \
    ${BASE_DIR}/aligned_Dfam${i}.txt \
    | sort -u -t '_' -k2 -n > ${BASE_DIR}/STAR_alignment_cons_dfam_${i}/AS.txt

    :> ${BASE_DIR}/aligned_ref_AS${i}.txt
    :> ${BASE_DIR}/aligned_inter${i}.txt
    :> ${BASE_DIR}/aligned_Rb${i}.txt

    echo "Annotate the components nodes with the TE and ref genes informations"
    echo "" > ${BASE_DIR}/aligned_Rb${i}.txt
    ##Annotate the components nodes with the TE and ref genes informations
    python3 ${WORK_DIR}/add_ref_TE.py \
            ${BASE_DIR}/comp${i}.txt \
            ${BASE_DIR}/aligned_ref_AS${i}.txt \
            ${BASE_DIR}/aligned_Dfam${i}.txt \
            ${BASE_DIR}/aligned_Rb${i}.txt \
            ${BASE_DIR}/aligned_inter${i}.txt \
            ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.abundance


    echo "Consensus sequence :" >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt

    python3 seq_consensium_of_comp.py \
        ${BASE_DIR}/comp${i}_annotated.nodes \
        ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
        ${K} >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt

    echo "Max abundance of the component :" >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt
    awk 'BEGIN {s=0} $9 > s {s=$9} END {print s}' FS="\t" ${BASE_DIR}/comp${i}_annotated.nodes >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt

    echo "Leaves abundance sum:" >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt
    awk '{s+=$9} END {print s}' FS="\t"  ${BASE_DIR}/genes_of_comp${i}/leaves.txt >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt


    echo "Analysis of the genes and TE"
    ##Analysis of the genes and TE
    ${WORK_DIR}/analysis_gene_TE \
            ${BASE_DIR}/comp${i}_annotated.nodes \
            ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges > \
            ${BASE_DIR}/gene_TE_analysis${i}.txt
done


pip install pandas
pip install matplotlib


mkdir -p ${BASE_DIR}/STAR_alignment_cons_dfam_all

#Convert the unitigs into reads for star
cp ${RESULTS_DIR}/graph/outputNodes.txt ${BASE_DIR}/all_unitigs.txt
python3 ${WORK_DIR}/reads_to_align.py ${BASE_DIR}/all_unitigs.txt ${BASE_DIR}/all_unitigs.fa 0

###--seedPerReadNmax 1000000 added because star asked to
${STAR_BIN_L} \
    --genomeDir ${BASE_DIR}/../dfam_ref \
    --runMode alignReads \
    --runThreadN 8 \
    --readFilesIn ${BASE_DIR}/all_unitigs.fa \
    --outFileNamePrefix ${BASE_DIR}/STAR_alignment_cons_dfam_all/ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFilterMultimapNmax 10000 \
    --winReadCoverageRelativeMin 0.5 \
    --seedMultimapNmax 10000000 \
    --alignWindowsPerReadNmax 100000 \
    --seedSearchStartLmax 1000 \
    --winAnchorMultimapNmax 50000 \
    --outReadsUnmapped Fastx


## STEP 2 : Analysis of the alignment and annotation of the nodes ###

echo "Filter the reads to only keep the best mapped ones"


python3 ${WORK_DIR}/filter_bam.py \
    ${BASE_DIR}/STAR_alignment_cons_dfam_all/Aligned.sortedByCoord.out.bam \
    ${BASE_DIR}/aligned_Dfam_all.txt

awk '{print $1"\t"$20}' FS=["\t:"] \
${BASE_DIR}/aligned_Dfam_all.txt \
| sort -u -t '_' -k2 -n > ${BASE_DIR}/STAR_alignment_cons_dfam_all/AS.txt


:> ${BASE_DIR}/aligned_Rb_all.txt
:> ${BASE_DIR}/aligned_ref_AS_all.txt
:> ${BASE_DIR}/aligned_inter_all.txt

echo "Annotate the components nodes with the TE and ref genes informations"
##Annotate the components nodes with the TE and ref genes informations
python3 ${WORK_DIR}/add_ref_TE.py \
        ${BASE_DIR}/all_unitigs.txt \
        ${BASE_DIR}/aligned_ref_AS_all.txt \
        ${BASE_DIR}/aligned_Dfam_all.txt \
        ${BASE_DIR}/aligned_Rb_all.txt \
        ${BASE_DIR}/aligned_inter_all.txt \
        ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.abundance


### STEP 3 : Analysis of the neighbours of the components ###

##The number of components to compute
MAXI=$(ls ${BASE_DIR}/comp*.txt | wc -l)
echo "Number of comps : ${MAXI}"

#touch ${BASE_DIR}/expressed_genes.nodes
#touch ${BASE_DIR}/expressed_genes.edges
for ((i=0; i<$MAXI; i++))
  do
  ## Dealing with the comp$i.txt file
  echo "Computing the genes of the component  ${i}..."
  mkdir -p ${BASE_DIR}/genes_of_comp$i

  ${WORK_DIR}/gene_finder.exe \
      ${BASE_DIR}/comp${i}.txt \
      ${BASE_DIR}/all_unitigs_annotated.nodes \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
      ${BASE_DIR}/genes_of_comp$i \
      ${K}

  echo "Consensus sequence :" >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt

  python3 seq_consensium_of_comp.py \
      ${BASE_DIR}/comp${i}_annotated.nodes \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
      ${K} >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt

  echo "Max abundance of the component :" >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt
  awk 'BEGIN {s=0} $9 > s {s=$9} END {print s}' FS="\t" ${BASE_DIR}/comp${i}_annotated.nodes >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt

  echo "Leaves abundance sum:" >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt
  awk '{s+=$9} END {print s}' FS="\t"  ${BASE_DIR}/genes_of_comp${i}/leaves.txt >> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt


  done

mkdir -p ${BASE_DIR}/genes_of_comp
./gene_finder.exe \
      ${BASE_DIR}/comp \
      ${BASE_DIR}/all_unitigs_annotated.nodes \
      ${BASE_DIR}/../graph/graph_hc_1_hc_2_k41_C0.05.edges \
      ${BASE_DIR}/genes_of_comp \
      ${K} \
      ${MAXI}

for ((i=0;i<${MAXI};i++)) do python3 analysis_comp.py ${BASE_DIR}/comp${i}.fa ${BASE_DIR}/aligned_ref_AS${i}.txt ${BASE_DIR}/aligned_Dfam${i}.txt ${BASE_DIR}/cr_ana${i}.txt --short; done
python3 CR_ana_comps.py ${BASE_DIR}/cr_ana ${MAXI} > ${BASE_DIR}/tot_cr_ana.txt
paste ${BASE_DIR}/genes_of_comp/transcript_summary_comps.txt ${BASE_DIR}/tot_cr_ana.txt > ${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated.txt



#Checking consensus sequences matching TE

: > ${BASE_DIR}/genes_of_comp/consensus_seqs.fa
for ((i=0; i<$MAXI; i++))
  do
  echo ">cons_${i}" >> ${BASE_DIR}/genes_of_comp/consensus_seqs.fa
  sed '8q;d' ${BASE_DIR}/genes_of_comp$i/gene_summary.txt >> ${BASE_DIR}/genes_of_comp/consensus_seqs.fa
  done

: > ${BASE_DIR}/genes_of_comp/consensus_seqs_ab.fa
for ((i=0; i<$MAXI; i++))
  do
  echo ">cons_${i}" >> ${BASE_DIR}/genes_of_comp/consensus_seqs_ab.fa
  awk '{print $9,$2}' FS="\t" ${BASE_DIR}/comp${i}_annotated.nodes | sort -k1,1nr | head -1 | awk '{print $2}' >> ${BASE_DIR}/genes_of_comp/consensus_seqs_ab.fa
  done


mkdir -p ${BASE_DIR}/STAR_alignment_seqs_consensus
${STAR_BIN_L} \
    --genomeDir ${BASE_DIR}/../dfam_ref \
    --runMode alignReads \
    --runThreadN 8 \
    --readFilesIn ${BASE_DIR}/genes_of_comp/consensus_seqs.fa \
    --outFileNamePrefix ${BASE_DIR}/STAR_alignment_seqs_consensus/ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFilterMultimapNmax 10000 \
    --winReadCoverageRelativeMin 0.5 \
    --seedMultimapNmax 10000000 \
    --alignWindowsPerReadNmax 100000 \
    --seedSearchStartLmax 1000 \
    --winAnchorMultimapNmax 50000 \
    --outReadsUnmapped Fastx

python3 ${WORK_DIR}/filter_bam.py \
    ${BASE_DIR}/STAR_alignment_seqs_consensus/Aligned.sortedByCoord.out.bam \
    ${BASE_DIR}/aligned_seqs_cons.txt

awk '{print $1"\t"$3"\t"$10}' FS=["\t"] \
${BASE_DIR}/aligned_seqs_cons.txt \
| sort -u -t '_' -k2 -n \
| cut -c6- > ${BASE_DIR}/STAR_alignment_seqs_consensus/seqs_cons_matching_TE.txt

awk '$2!="*"{print $1}' FS=["\t"] ${BASE_DIR}/STAR_alignment_seqs_consensus/seqs_cons_matching_TE.txt > ${BASE_DIR}/id_comp_matching_TE.txt


awk 'NR==FNR {ids[$1]; next}
{
    idx = FNR - 1
    if ($5=="TE" && (idx in ids))        print $0 "\tTP";
    else if ($5=="TE" && !(idx in ids))  print $0 "\tFN";
    else if ($5!="TE" && (idx in ids))   print $0 "\tFP";
    else                                 print $0 "\tTN";
}' FS=["\t"] ${BASE_DIR}/id_comp_matching_TE.txt \
${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated.txt > ${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated_cons.txt

awk '$6=="TP"{pos+=1} $6=="FN"{faux+=1} END {print pos, faux}' FS=["\t"]  \
${BASE_DIR}/genes_of_comp/transcript_summary_comps_annotated_cons.txt > ${BASE_DIR}/genes_of_comp/seq_cons_precision.txt


## TeTools to know the count of each TE
awk '$1 ~ /^>/ {print $0, $2}' FS=['\t'_] ${DFAM_FA} | cut -c 2- > ${SPE_DIR}/rosette.fa



${TECOUNT} -rosette ${SPE_DIR}/rosette.fa  \
    -column 2 \
    -TE_fasta ${DFAM_FA} \
    -count ${SPE_DIR}/count_TE.txt \
    -bowtie2 \
    -RNA ${READS_1} \
    -RNApair ${READS_2} \
    -insert 500