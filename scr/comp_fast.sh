#!/bin/bash
#Need some comments here


begin=`date +%s`

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

if [[ ${SKIP_BUILD_RUST} == "" ]]; then
  ##Build the bin for the graph
  echo "Building the bin for the graph..."
  cargo build --release --manifest-path ${WORK_DIR}/graph/Cargo.toml
  cp ${WORK_DIR}/graph/target/release/graph ${WORK_DIR}/graph.exe

  ##Build the bin for the gene_finder function
  echo "Building the bin for the gene_finder function..."
  cargo build --release --manifest-path ${WORK_DIR}/gene_finder_de_novo/Cargo.toml
  cp ${WORK_DIR}/gene_finder_de_novo/target/release/gene_finder_de_novo ${WORK_DIR}/gene_finder_de_novo.exe


  ##Build the bin for the filtering_low_ab_percent function
  echo "Building the bin for the filtering_low_ab_percent function..."
  cargo build --release --manifest-path ${WORK_DIR}/filtering_low_ab_percent/Cargo.toml
  cp ${WORK_DIR}/filtering_low_ab_percent/target/release/filtering_low_ab_percent ${WORK_DIR}/filtering_low_ab_percent.exe

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Build time (in seconds): $elapsed \n"
fi


if [[ ${SKIP_FASTP_HC} == "" ]]; then
  ##FastP of the reads to remove the poly(A) tails
  echo "FastP of the reads ..."

  ${FASTP_BIN} \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --trim_poly_x \
      --thread ${P} \
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

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "FastP time and HC (in seconds): $elapsed \n"
fi


if [[ ${SKIP_GEN_GRAPH} == "" ]]; then

  ##Compute the DGB with bcalm
  echo "DGB with bcalm ..."
  ${BCALM_BIN} \
      -in ${DATA_DIR}/${SPE}/hc_1.fq \
      -in ${DATA_DIR}/${SPE}/hc_2.fq \
      -kmer-size ${K} \
      -abundance-min 2 \
      -nb-cores ${P} \
      -out-dir ${RESULTS_DIR}/graph


  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "DGB time (in seconds): $elapsed \n"


  ##Filter the low abundance percent unitigs
  echo "Filtering the low abundance percent unitigs ..."
  ${WORK_DIR}/filtering_low_ab_percent.exe \
      ${RESULTS_DIR}/graph/hc_1_hc_2_k${K}.nodes \
      ${RESULTS_DIR}/graph/hc_1_hc_2_k${K}.edges \
      ${RESULTS_DIR}/graph/hc_1_hc_2_k${K}.abundance \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K} \
      5

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Filtering time (in seconds): $elapsed \n"


  ##Compute the weighting
  echo "Weighting of the nodes..."
  g++ -g graph.cpp ponderation.cpp -o graph.exe
  ${WORK_DIR}/graph.exe \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.nodes \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
      ${D_NT} \
      -k ${K}  \
      -o ${RESULTS_DIR}/graph/outputNodes.txt


  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Weighting time (in seconds): $elapsed \n"


  ##Compute the threshold
  echo "Threshold of the nodes..."
  T=$(python3 ${WORK_DIR}/plot.py ${RESULTS_DIR}/graph/outputNodes.txt top1)
  # Update the T variable in the environment.sh file
  sed -i "s/^T=.*/T=${T}/" "${ENV}"
  echo "T=${T}"

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Threshold time (in seconds): $elapsed \n"

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

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Agglomeration time (in seconds): $elapsed \n"

fi


##The number of components to compute
MAXI=$(ls ${BASE_DIR}/comp*.txt | wc -l)
echo "Number of comps : ${MAXI}"

##Activation the virtual environment for python
source ${WORK_DIR}/venv/bin/activate

##Compute the analysis over every component
echo "Consensus sequence :" #>> ${BASE_DIR}/genes_of_comp${i}/gene_summary.txt
python3 seq_consensium_of_comps.py \
    ${BASE_DIR}/comp \
    ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
    ${BASE_DIR}/seq_consensium.txt \
    ${K} \
    ${MAXI}

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Consensus sequences time (in seconds): $elapsed \n"

${WORK_DIR}/gene_finder_de_novo.exe \
      ${BASE_DIR}/comp \
      ${RESULTS_DIR}/graph/outputNodes.txt \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}_C0.05.edges \
      ${RESULTS_DIR}/graph/graph_hc_1_hc_2_k${K}.abundance \
      ${BASE_DIR}/genes_of_comp \
      ${K} \
      ${MAXI}

python3 analysis_comp_de_novo.py \
    ${BASE_DIR}/comp \
    ${BASE_DIR}/seq_consensium.txt \
    ${BASE_DIR}/analysis_comp

  end=`date +%s`
  elapsed=`expr $end - $begin`
  begin=`date +%s`
  echo "Gene finding time (in seconds): $elapsed \n"
