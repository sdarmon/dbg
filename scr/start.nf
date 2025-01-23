#!/usr/bin/env nextflow

workDir = ''

println "\tDébut de la Pipeline Nextflow\nA executer dans le dossier scr.\nExecution type: \n\t-nextflow run tutorial.nf --path [YourWorkDirPath]\n\t-pwd | xargs nextflow run tutorial.nf --path : Permet de récupérer directement le WorkDir dedans\n "


// Pour le bon fonctionnement de ce programme, il est important d'indiquer les fichiers par leur chemins absolus.
// La variable workDir contient par défaut le chemin du dossier de travail.

if (params.path != null){
workDir = params.path
println "Path loaded\n"
}

// Pour chaque jeu de données à traiter, il faut reprendre le modèle ci-dessous, permettant d'indiquer les chemins
// des différents fichiers.

test = ["name":"", "nodes":"","abundance":"",  "edges":"", "rayon":"", "k" : "", "threshold" : ""]
test["name"] = "test"
test["nodes"] = "${workDir}/../data/test.nodes"
test["edges"] = "${workDir}/../data/test.edges"
test["abundance"] = "${workDir}/../data/test.abundance"
test["rayon"] = 3   //Pour vrai jeu de donnée prendre 10
test["k"] = 5   //Pour vrai jeu de donnée prendre 41
test["threshold"] = 4   //Pour vrai jeu de donnée prendre 14

// Enfin, pour lancer l'execution des différents jeux de données, il suffit de mettre leur noms à l'intérieur Channel suivant :

a_executer = Channel.from(test)  // Marche avec plusieurs fichiers, avec 3 cela donnerait : "Channel.from(test1,test2,test3)"

// Il n'y a plus rien à modifier, il n'y a plus qu'à lancer le NewtFlow
// à l'aide de la commande `pwd | xargs nextflow run tutorial.nf --path `

process calculpoids {
    input:
    val spe from a_executer

    output:
    val spe into pondere

    script:
    name = spe.name
    nodes = spe.nodes
    edges = spe.edges
    r = spe.rayon
    k = spe.k
    c = spe.threshold
    """
    g++ ${workDir}/graph.cpp ${workDir}/ponderation.cpp -o ${workDir}/ponderation.exe
    ${workDir}/ponderation.exe  ${nodes} ${edges} ${r} -k ${k} -o ${workDir}/../processing/${name}.nodes.pondere
    python3 ${workDir}/reads_to_align.py ${workDir}/../processing/${name}.nodes.pondere ${workDir}/../processing/${name}_clean.nodes.pondere 0 -clean
    """
}


process gene_composantes {
    input:
    val spe from pondere

    output:
    val spe into simpl

    script:
    name = spe.name
    nodes = spe.nodes
    edges = spe.edges
    r = spe.rayon
    k = spe.k
    c = spe.threshold
    """
    g++ ${workDir}/graph.cpp ${workDir}/gene_comp.cpp -o ${workDir}/gene_comp.exe
    ${workDir}/gene_comp.exe ${workDir}/../processing/${name}_clean.nodes.pondere ${edges} -c ${c} -k ${k} ${workDir}/..
    """
}


process simplification {
    input:
    val spe from simpl


    script:
    name = spe.name
    nodes = spe.nodes
    edges = spe.edges
    r = spe.rayon
    k = spe.k
    c = spe.threshold
    """
    g++  ${workDir}/graph.cpp ${workDir}/composantes.cpp -o ${workDir}/composantes.exe
    ${workDir}/composantes.exe ${nodes} ${edges} -k ${k} ${workDir}/../processing/comp ${workDir}/../results/clean
    """
}