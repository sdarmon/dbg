#!/bin/bash
#

# Fonction pour afficher l'aide
afficher_aide() {
    echo "Usage: $0 [-h] graph.nodes graph.edges graph.ab comp.txt out.txt d ab_min"
    echo
    echo "This function extend a comp up to its d-th neighbors while annotating all the nodes."
    echo
    echo "Options:"
    echo "  -h    Show this command"
    echo
    echo "Arguments:"
    echo "  graph.nodes        The node file of the whole graph"
    echo "  graph.edges        The edge file of the whole graph"
    echo "  graph.ab           The abundance file of the whole graph"
    echo "  comp.txt           The file of the comp"
    echo "  all_unitigs_annotated.nodes            The file of the all unitigs annotated"
    echo "  dir_out            Output file directory"
    echo "  d                  distance around the comp"
    echo "  ab_min             Minimum abundance of the nodes to be kept"

    echo
    echo "Exemple:"
    echo "  $0 \$\{BASE_DIR\}/../../mou/graph/graph_hc_1_hc_2_k41.nodes \$\{BASE_DIR\}/../../mou/graph/graph_hc_1_hc_2_k41_C0.05.edges \$\{BASE_DIR\}/../../mou/graph/graph_hc_1_hc_2_k41.abundance \$\{BASE_DIR\}/comp3_annotated.nodes \$\{BASE_DIR\}/all_unitigs_annotated.nodes ~/Documents/data/graph/comp3 10"
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
            echo "Option invalide : -$OPTARG" >&7
            afficher_aide
            ;;
    esac
done

# Supprimer les options traitées des arguments
shift $((OPTIND - 1))

# Vérification des arguments restants
if [ "$#" -lt 7 ]; then
    echo "Erreur : Nombre insuffisant d'arguments."
    afficher_aide
fi

# Variables
nodes=$1
edges=$2
ab=$3
in=$4
all_unitigs=$5
out=$6
d=$7
base=$(basename "$out")
#Computing the neighborhood

./neigh.exe ${nodes} ${edges} ${in} -o ${out}_n -c 0 -d "$d"

### Get the nodes IDs
awk '{print $1}' FS="\t" ${out}_n.nodes | sort -u -n > ${out}_n_IDs.txt

### Get the lines from the in file such that the node ID is in the list of IDs
awk 'NR==FNR {ids[$1]; next} $1 in ids {print}' FS="\t" ${out}_n_IDs.txt ${all_unitigs} > ${out}_n_annotated.nodes

### Sort the edges file to suppress the duplicates
sort -u ${out}_n.edges > ${out}_n_sorted.edges

