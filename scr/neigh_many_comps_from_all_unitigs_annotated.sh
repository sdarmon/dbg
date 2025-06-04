#!/bin/bash
#/home/sdarmon/pbil_home/Documents/dbg/scr/neigh_many_comps_from_all_unitigs_annotated.sh
# Fonction pour afficher l'aide
afficher_aide() {
    echo "Usage: $0 [-h] graph.nodes graph.edges graph.ab prefix_comp prefix_out d nb_comp ab_min"
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
    echo "  comp.txt           prefix of the comps"
    echo "  all_unitigs_annotated.nodes            The file of the all unitigs annotated"
    echo "  dir_out            Output file directory"
    echo "  d                  distance around the comp"
    echo "  nb_comp            Number of components to be extended"
    echo "  ab_min             Minimum abundance of the nodes to be kept"

    echo
    echo "Exemple:"
    echo "  $0 \$\{BASE_DIR\}/../../mou/graph/graph_hc_1_hc_2_k41.nodes \$\{BASE_DIR\}/../../mou/graph/graph_hc_1_hc_2_k41_C0.05.edges \$\{BASE_DIR\}/../../mou/graph/graph_hc_1_hc_2_k41.abundance \$\{BASE_DIR\}/comp \$\{BASE_DIR\}/all_unitigs_annotated.nodes ~/Documents/data/graph/comp_all 100 100"
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
            echo "Option invalide : -$OPTARG" >&9
            afficher_aide
            ;;
    esac
done

# Supprimer les options traitées des arguments
shift $((OPTIND - 1))

# Vérification des arguments restants
if [ "$#" -lt 9 ]; then
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
n=$8
ab_min=$9
base=$(basename "$out")
#Computing the neighborhood

#Looping over the n components from 0 to n-1
for ((i=0; i<$n; i++))
 do
    echo "Computing the neighborhood of component $i"
    ./neigh.exe ${nodes} ${edges} ${in}${i}_annotated.nodes -o ${out}_${i} -c 0 -d "$d"
done


### Get the nodes IDs for the components #awk -v ab="${ab_min}" ' $9 >= ab {print $1}' FS="\t"
echo "" > ${out}_IDs.txt
echo "" > ${out}_IDs_comps.txt
for ((i=0; i<$n; i++))
     do
    awk '{print $1}' FS="\t" ${out}_${i}.nodes >> ${out}_IDs.txt
    awk '{print $1}' FS="\t" ${in}${i}_annotated.nodes  >> ${out}_IDs_comps.txt
done

sort -u ${out}_IDs.txt > ${out}_IDs_sorted.txt

### Get the lines from the in file such that the node ID is in the list of IDs
awk -v ab="${ab_min}"  'NR==FNR {ids[$1]; next} $1 in ids && $9 >= ab {print}' FS="\t" ${out}_IDs_sorted.txt ${all_unitigs} > ${out}_n.nodes
awk 'NR==FNR {ids[$1]; next} $1 in ids {print}' FS="\t" ${out}_IDs_comps.txt ${all_unitigs} >> ${out}_n.nodes

sort -u -k1,1n ${out}_n.nodes > ${out}_annotated.nodes
### Merge the edges files and sort to suppress the duplicates

awk 'NR==FNR {ids[$1]; next} $1 in ids && $2 in ids && $1 < $2 {print}' FS="\t" ${out}_annotated.nodes ${edges} > ${out}_n.edges
sort -u ${out}_n.edges > ${out}_sorted.edges

#Remove all the temporary files
echo "Removing temporary files"
rm ${out}_IDs.txt ${out}_IDs_sorted.txt ${out}_IDs_comps.txt
rm ${out}_n.nodes
rm ${out}_n.edges
for ((i=0; i<$n; i++))
    do
    rm ${out}_${i}.nodes ${out}_${i}.edges
    done
