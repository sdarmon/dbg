#!/bin/bash
#Need some comments here



#Loading the environment variables and paths through
#the first argument or the default environment.sh
if [[ -n "$1" && "$1" == *.sh ]]; then
    source "$1"
else
    source environment.sh
fi

#Create the general directories of that species
mkdir -p $DATA_DIR
mkdir -p $RESULTS_DIR
mkdir -p $SPE_DIR
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
    --correction \
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
python3 homomorphic_compression.py  ${READS_1}.fastp ${DATA_DIR}/hc_1.fq 5
python3 homomorphic_compression.py  ${READS_2}.fastp ${DATA_DIR}/hc_2.fq 5

##Compute the TEs

# Extraction des TEs de la base de données Dfam
echo "Extraction of the TEs from Dfam ..."
source ${VENV_FAMDB}
${FAMDB_BIN} -i ${LIBRARY_DIR} families \
--include-class-in-name \
--curated \
--ancestors \
--descendants \
"${SPE_NAME}" --format fasta_name \
 | sed 's/#/\t/g' \
 | sed 's/ @/\t/g' \
 | sed 's/^>/>dfam_/g' > ${DFAM_FA}


# Extraction des TEs de la base de données Repbase
echo "Extraction of the TEs from Repbase ..."
${FAMDB_BIN} -i ${LIBRARY_DIR} lineage -ad "${SPE_NAME}" \
| awk -F '[0-9<(]+' 'NR>2 {print $2}' \
| sed 's/^ //g' \
| sort -u > ancestors.txt

awk 'NR==FNR {ancestors[$1]; next} /^>/ {split($0, fields, "\t"); species=fields[3]; print_seq = (species in ancestors)} print_seq' ancestors.txt ${REPEATMASKER_DIR}/repBase2508.fasta \
| sed 's/^>/>rb_/g' \
| awk '{if($0 ~ /^>/) {print $0} else {print toupper($0)}}'  > ${RB_FA}

#Nettoyage
rm ancestors.txt



##Compute the DGB with kissplice
ulimit -s unlimited
echo "DGB with kissplice ..."
$KISSPLICE_BIN \
    -r ${DATA_DIR}/hc_1.fq \
    -r ${DATA_DIR}/hc_2.fq \
    -k ${K} \
    -o ${RESULTS_DIR}/graph

##Compute the weighting
echo "Weighting of the nodes..."
# g++ -g graph.cpp ponderation.cpp -o graph.exe
${WORK_DIR}/graph.exe \
    ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.nodes \
    ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
    ${D_NT} \
    -k ${K}  \
    -o ${RESULTS_DIR}/graph/outputNodes.txt

## Compute the connexe components
echo "Agglomeration of connexe components..."
g++ -g ${WORK_DIR}/graph.cpp ${WORK_DIR}/agglo.cpp -o ${WORK_DIR}/agglo.exe

rm ${BASE_DIR}/*
${WORK_DIR}/agglo.exe \
    ${RESULTS_DIR}/graph/outputNodes.txt \
    ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
    -c ${T} \
    -d ${D_NT} \
    ${RESULTS_DIR} \
    -clean ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.abundance > ${RESULTS_DIR}/rapportAgglo.txt
#agglo.exe ${RESULTS_DIR}/graph/outputNodes.txt  ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges -c ${T} -d ${D_NT} ${RESULTS_DIR} -clean ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.abundance > ${RESULTS_DIR}/rapportAgglo.txt

## Compute the genome STAR ref
echo "Generating genome STAR ref..."

$STAR_BIN \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir $RESULTS_DIR/ref \
--genomeFastaFiles $GEN_FA \
--sjdbGTFfile $GEN_GTF \
--genomeSAindexNbases 13 \
--genomeSAsparseD 4 \
--genomeChrBinNbits 16

$STAR_BIN \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ${RESULTS_DIR}/dfam_ref \
--genomeFastaFiles $DFAM_FA \
--genomeSAindexNbases 9 \
--genomeSAsparseD 4 \
--genomeChrBinNbits 16

$STAR_BIN \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ${RESULTS_DIR}/rb_ref \
--genomeFastaFiles $RB_FA \
--genomeSAindexNbases 9 \
--genomeSAsparseD 4 \
--genomeChrBinNbits 16

##Build the bin for gene and TE analysis
cargo build --release --manifest-path ${WORK_DIR}/graph/Cargo.toml
cp ${WORK_DIR}/graph/target/release/graph ${WORK_DIR}/analysis_gene_TE

##Build the bin gene_finder function
cargo build --release --manifest-path ${WORK_DIR}/gene_finder/Cargo.toml
cp ${WORK_DIR}/gene_finder/target/release/gene_finder ${WORK_DIR}/gene_finder.exe
