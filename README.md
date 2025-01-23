# Stage M2 : Développements de méthodes d’analyse de l’épissage alternatif à partir de données RNA-seq

## Présentation

Ce projet contient contient les résultats de mon stage de M2 portant sur le développement de méthodes d’analyse de l’épissage alternatif à partir de données RNA-seq. Un rapport de stage est disponible sur ce git. L'algorithme principal de ce Git permet de simplifier les graphes de De Bruijn afin de trouver de nouveaux évènements d'épissage alternatif.

Voici une présentation du contenu de ce Git :

### Data :

Ce dossier contient un mini jeu de données pour tester les algorithmes mis en place. Ce jeu de donnée "test" correspond à un graphe de De Bruijn compacté et est décomposé en trois fichiers :

- `test.nodes` : Chaque ligne correspond à un unique sommet qui est défini par un numéro et une séquence de nucléotide.
- `test.edges` : Chaque ligne correspond à une unique arête qui est définie par un sommet de départ, un sommet d'arrivé et un type d'arête (plus d'explications sur les types des arêtes à la structure du graphe sont disponibles dans le rapport).
- `test.abundance` : Chaque ligne correspond à l'abondance des k-mers de chaque sommet. Je n'utilise pas cette métrique, elle ne sert qu'au bon fonctionnement de *KisSplice*.

### Processing :

- `test.nodes.pondere` : Ce fichier correspond au fichier `test.nodes` mais avec en plus pour chaque sommet son poids. Il a était obtenu à la suite de l'étape 1 de l'algorithme.
- `comp0.txt` : Ce fichier contient tous les sommets faisant parties de la composante 0. Il a été généré à la suite de l'étape 2.

### Results :

Ce dossier contient les résultats finaux de mon alogithme. Il est composé de fichiers `clean.nodes`, `clean.edges` et `clean.abundance` qui correspondent au nouveau graphe qui sera par la suite utilisé par *KisSplice*.


### Scr :

Ce dossier contient tous les algorithmes. On peut les décomposer en trois catégories :

Tout d'abord, il y a les programmes de l'algorithme de simplification : 

- `graph.cpp` et `graph.h` : Contiennent toutes les structures et fonctions liées aux graphes que j'ai utilisés.
- `ponderation.cpp` : Etape 1 de l'algorithme permettant de faire la pondération du graphe.
- `gene_comp.cpp` : Etape 2 de l'algorithme permettant de générer toutes les composantes.
- `composantes.cpp` : Etape 3 de l'algorithme permettant de simplifier effectivement le graphe.

Ensuite, il y a les programmes liés à l'analyse de mes données ou de mes méthodes :

- `empty_composantes.cpp` : Ce programme fonctionne comme `composantes.cpp`, supprime les composantes, mais ne rajoute pas les nouvelles arêtes! Cela m'a permis d'avoir des graphes de contrôle.
- `filtrage_bulle.py` : Cette fonction permet d'annoter les bulles afin de savoir les quelles sont nouvelles ainsi que de savoir quelles bulles passent par quelles composantes.
- `interBulle.py` : Cette fonction permet de garder qu'une bulle par composante, et est utilisée afin d'avoir un aperçu des composantes.
- `neighborhood.cpp` : Ce programme permet d'extraire un sous-graphe à distance donnée, permettant d'étudier des sous-graphes.
- `plot.py` : Cette fonction me permet de produire divers graphiques (dont notamment ceux de mon rapport) facilitant l'analyse.
- `rapportAgglo.py` : Cette fonctionne permet d'intersecter les comosantes avec des éléments transposables connues et d'afficher un graphique récapitulatif du contenu des composantes. De plus, elle affiche également comment sont répartis tous les éléments transposables connus (en les classant en fonction de dans combien de composantes ils sont). 
- `reads_to_align.py` : Cette fonction permet de préparer le fichiers des sommets du graphes de plusieurs manières; soit elle va nettoyer les données (queue poly(A)) ou alors elle peut transformer la liste des sommets en fichier `.fa`, alignable sur des génomes de référence. 
- `suppDoublon.py` : Cette fonction permet de filtrer des éléments transposables en double ou des bulles en double.

Finalement, il y a les programmes de pipeline utilisant le logiciel *Nextflow*:

- `start.nf` : Permet d'exécuter l'algorithme sur n'importe quel graphe en un ligne. Par défault, s'exécute uniquement sur le graphe `test`.
- `analyse_TE` : Permet d'effectuer l'analyse complète des éléments transposables des jeux de données. Nécéssite d'avoir des fichiers d'annotation pour fonctionner.




## Exemple d'utilisation 

A partir du logiciel *Nextflow*, il est possible d'exécuter toutes les étapes de l'algorithme sur n'importe quel graphe. Par défaut, `start.nf` s'exécute sur le graphe `test`, mais il est possible de rajouter n'importe quel graphe à traiter directement dans le fichier `start.nf`.

Pour exécuter le programme, il suffit d'exécuter la commande suivante :

```
pwd | xargs nextflow run start.nf --path
``` 

Sinon, il est également possible de refaire chaque étape à la main :

### Etape 1 : la podération du graphe

```
g++ graph.cpp ponderation.cpp -o ponderation.exe 
./ponderation.exe ../data/test.nodes ../data/test.edges 3 -k 5 -o ../processing/test.nodes.pondere
python3 reads_to_align.py ../processing/test.nodes.pondere ../processing/test_clean.nodes.pondere 0 -clean
``` 

### Etape 2 : la génération des composantes

```
g++ graph.cpp gene_comp.cpp -o gene_comp.exe 
./gene_comp.exe ../processing/test_clean.nodes.pondere ../data/test.edges -c 4 -k 5 ..
``` 

### Etape 3 : la simplification du graphe de De Bruijn

```
g++ graph.cpp composantes.cpp -o composantes.exe 
./composantes.exe ../processing/test_clean.nodes.pondere ../data/test.edges -k 5 ../processing/comp0.txt 1 ../results/clean
``` 

Ainsi, on obtient le graphe simplifié dans le dossier `results` qui peut être ensuite utilisé par le logiciel *KisSplice* pour obtenir une énumération de bulles.




## Améliorations à apporter

### Sur l'étape 2:

J'effectue un tri des sommets car historiquement, je ne conservais que les 100 composantes ayant les poids les plus élevés. Cependant, avoir ce tri sur les composantes n'est pas forcément pertinant, et n'est utile que pour l'analyse des composantes. Il serait donc préférable d'utiliser l'algorithme de Tarjan pour obtenir les composantes linéairement lors d'une implémentation concrète.


### Sur l'étape 3:

J'ai décidé de choisir pour représenter l'arête consensus, de prendre le plus court chemin entre les sommets à relier. Ce choix là n'est pas forcément utile et le fait de relacher cette condition permet de concevoir un algorithme linéaire en les tailles des données. Etant actuellement l'étape limitante, il s'agit là de la meilleure optimisation à faire.


### Sur *KisSplice*

Pendant l'énumération des bulles, il peut être intéressant de prendre en compte les composantes. En effet, on ne souhaite pas énumérer les bulles dont les deux chemins passent par la même composante. Cela permetterait alors d'améliorer l'efficacité de l'énumération.
