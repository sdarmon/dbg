# Scripts Utiles

Voici la liste des différentes commandes que j'ai effectué afin d'obtenir différentes analyses.
Attention, les noms ou les chemins des fichiers peuvent avoir changé. De plus, comme je faisais beaucoup de modifications au fur et à mesure, tous les fichiers ne sont pas forcément compatibles; et dans le doute, il vaut mieux tous relancer. Pour toutes les fonctions, les exécuter sans argument affichera les arguments classiques à ajouter.

## Composition des dossiers

### data

Ce dossier est divisé en sous dossier pour chacun des jeux de données et contient les données principales.

### results

Ce dossier contient tous les résultats obtenues par les différentes méthodes et algorithmes.
Le dossier le plus important est le dossier `moustique` qui est divisé en plusieurs types de dossiers :

- `compX` : ces dossiers contiennent les composantes obtenues à la fin de l'algo `gene_comp.exe`. J'ai supprimé sans faire exprès la comp10, il faut la reformer en cas d'utilisation.
- `graphCompX` : ces dossiers contiennent le graphe final (à la suite de `composantes.exe`) et tous les autres fichiers finaux.
- `emptyGraph10BX` : ces dossiers contiennent le graphe final sans les composantes (à la suite de `empty_composantes.exe`) et tous les autres fichiers finaux, avec le paramètre `B` étant le branchement de `KisSplice`.

### stage-M2

Il s'agit d'un copie du Git, il y a tous les scripts à l'intérieur.

## Pour obtenir le nombre de composante : 

```
ls -l comp* | wc -l
```

### Exemple de la dernière étape de simplification du graphe qu'il faudrait relancer

```
g++ graph.cpp composantes.cpp -o composantes.exe
./composantes.exe ../../data/moustique/outputGraphmoustiqueClean.txt ../../../kissplice_results/kissplice_moustique/graph_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_C0.05.edges -k 41 ../../results/moustique/comp14/processing/comp 42806 ../../results/moustique/graphComp14/clean
```

## Analyse du moustique 

Alignement des bulles : 
```
STARlong --genomeDir ../DmGoth/results/moustique/genome/ --readFilesIn ../DmGoth/results/moustique/graphComp14/results_clean_type_1.fa
```

Kissplice2refgenome :
```
./kissplice2refgenome -a ../peda/DmGoth/data/moustique/Aedes_aegypti_lvpagwg.AaegL5.53.gff3 --count 0 ../peda/kissplice2refgenome/Aligned.out.sam
```

On génère le fichier `interBulle.txt` qui contient l'annotation des bulles et compare avec les bulles de références
```
python3 interBulle.py ../../results/moustique/comp14/processing/ ../../results/moustique/graphComp14/results_clean_type_1.fa 42606 14 -k 41 -compare ../../../kissplice_results/kissplice_moustique/results_IR03_B_R1_IR03_C_R1_IR03_D_R1_IR03_E_R1_IR13_B_R1_IR13_C_R1_IR13_D_R1_IR13_E_R1_k41_coherents_type_1.fa > ../../results/moustique/graphComp14/interBulle.txt
```

On génère `bulleUniqCompVrai.fa` qui contient une bulle par composante
```
python3 filtrageBulle.py ../../results/moustique/comp14/processing/  ../../results/moustique/graphComp14/results_clean_type_1.fa 42606 14 -k 41 -rapport > ../../results/moustique/graphComp14/bulleUniqComp.fa
```

On enlève la ligne d'annotation :
```
python3 quatreSurCinq.py ../../results/moustique/graphComp14/bulleUniqComp.fa > ../../results/moustique/graphComp14/bulleUniqCompVrai.fa
```

On re-aligne avec STAR puis analyse avec kiss3refgenome :
```
STARlong --genomeDir ../DmGoth/results/moustique/genome/ --readFilesIn ../DmGoth/results/moustique/graphComp14/bulleUniqCompVrai.fa

./kissplice2refgenome -a ../peda/DmGoth/data/moustique/Aedes_aegypti_lvpagwg.AaegL5.53.gff3 --count 0 ../peda/kissplice2refgenome/Aligned.out.sam
```

Ensuite, en filtrant les codes CIGAR de `Aligned.out.sam` on peut identifier des insertions d'éléments transposables ou on peut directement regarder le fichier `event.tsv` de kissplice2refgenome.




## Un autre moyen pour comparer les évènements d'épissage alternatif et en ragardant les chemins du bas :

```
python3 under_path.py ../../results/moustique/graphComp14/results_clean_type_1.fa ../../results/moustique/graphComp14/underPathCounted.txt
awk {'print $1'} ../../results/moustique/graphComp14/underPathCounted.txt | sort > ../../results/moustique/graphComp14/underPathSorted.txt
```

En faisant ça sur différents fichiers de bulles, on peut générer pleins de `underPathSorted` puis faire des `diff` dessus afin d'obtenir les différences.


## Obtenir uniquement les nouvelles bulles :

 Par exemple, pour obtenir les nouvelles bulles, on peut faire :

```
grep "only in" ../../results/moustique/graphComp14/interBulle.txt | awk {'print $(NF-3)'} > ../../results/moustique/graphComp14/id_new.txt
cat ../../results/moustique/graphComp14/id_new.txt | while read p; do grep "${p}" ../../results/moustique/graphComp14/results_clean_type_1.fa -A 1 > new_bulles.fa
```

Et on peut ensuite répéter le point précédent sur `new_bulles.fa` pour avoir ces chemins du bas.


## Intersection et trie des alignements de TE avec les alignements des unitigs :
 
```
bedtools intersect -split -a ../../data/droso/TE.gtf -b ../../results/droso/STAR_alignment/Aligned.sortedByCoord.out.bam -wo > ../../results/droso/intersectionBloc.txt
sort -k19 -n -r ../../results/droso/intersectionBloc.txt > ../../results/droso/intersectionSorted.txt```
```


## En cas de soucis :

Ne pas hésiter à m'envoyer un mail!