/*
 * Ce programme permet de faire une agglomeration des composantes
 * des graphes.
 */

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <map>
#include <queue>
#include "graph.h"
#define MAX 1024


//Fonction permettant d'enregistrer une composante du graphe `G` dans un fichier
void save_comp(Graph &G, vector<int> &compo, string outputPrefix, int rang){
    ofstream output;
    output.open(outputPrefix+"/processing/comp"+to_string(rang)+".txt");
    for (vector<int>::iterator it = compo.begin(); it != compo.end(); it++){
        output << *it << "\t" << G.Vertices[*it].label << "\t" << G.Vertices[*it].weight  << "\n";
    }
    return;
}


struct indexDic {
    int key;
    int weight;
    friend bool operator< (indexDic const& lhs, indexDic const& rhs) {
        return (lhs.weight < rhs.weight);
    }
};

int main(int argc, char** argv) {
    if (argc != 8) {
        cout << "Expected use of this program: \n\n\t" << argv[0]
             << " file.nodes file.edges -c value -k kmer outputPrefix \n" << endl;
        return 0;
    }

    //On charge le graphe
    vector <Edge> E;
    vector <Node> V;

    string nodesPath = argv[1];
    string edgesPath = argv[2];
    string outputPrefix;
    ifstream edges(edgesPath, std::ios::binary);
    ifstream nodes(nodesPath, std::ios::binary);

    read_edge_file(edges, E);
    read_node_file_weighted(nodes, V);

    Graph G(V, E);

    G.kmer = stoi(argv[6]);
    outputPrefix = argv[7];

    cout << "Graphe chargé et construit" << endl;

    int index;
    //Construction du vector `vu_total`
    vector<int> vu_total(G.N, 0);

    vector <vector<int>> components;
    components.clear();

    //Ici le critère choisi pour définir une composante est juste d'être le sous-graphe connexe maximum du graphe G
    //auquel on a enlevé tous les sommets de poids inférieur au threshold `threshold`.
    int threshold;
    threshold = atoi(argv[4]);
    int m = 0;
    int compt;
    vector<int> vu(G.N, 0);
    set<int> setVu;
    queue < Neighbor * > aVoir;

    cout << "Début du tas min" << endl;
    priority_queue <indexDic> indexTrie;
    for (int sommet = 0; sommet < G.N; sommet++) {
        if (G.Vertices[sommet].weight >= threshold) {
            indexTrie.push(indexDic{sommet, G.Vertices[sommet].weight});
        }
    } //Attention ici il y a une opti possible : on est en N log N mais c'est clairement possible de faire en N
    index = indexTrie.top().key;
    indexTrie.pop();
    vu_total[index] = 1;

    cout << "Début de la recherche des composantes" << endl;
    while (!indexTrie.empty()) //On cherche des composantes tant qu'il existe encore un sommet vérifiant
        //le critère
    {
        //On réalise alors un BFS depuis le sommet `index`
        for (set<int>::iterator it = setVu.begin(); it != setVu.end(); it++) {
            vu[(*it)] = 0;
        }
        setVu.clear();
        vector<int> compo;
        compo.clear();
        //aVoir étant une queue de BFS; est toujours censé être vide ici. Pour la sanité des algos qui suivent on la
        //vide au cas où...
        while (!aVoir.empty()) {
            aVoir.pop();
        }
        vu[index] = 1;
        setVu.insert(index);
        for (vector<Neighbor>::iterator it = G.Neighbors(index)->begin(); it != G.Neighbors(index)->end(); ++it) {
            if (G.Vertices[it->val].weight >= threshold) {
                aVoir.push(&(*it));
            }
        }
        G.BFS_func(threshold, aVoir, vu, setVu);

        //Puis on ajoute les sommets trouvés à la composante
        compt = 1;
        compo.push_back(index);
        vu_total[index] = 1;
        for (int i = 0; i < G.N; i++) {
            if (vu[i] and i != index) {
                compt++;
                compo.push_back(i);
                vu_total[i] = 1;
            }
        }
        if (compt > 1) {
            save_comp(G, compo, outputPrefix, m); //On enregistre la composante sur l'ordinateur
            m++;
            components.push_back(compo); //Et on ajoute la composante au vecteur de composantes
            cout << "Composante trouvée de départ " << index << " et de poids " << G.Vertices[index].weight
                 << " et de taille " << compo.size() << endl;
        }

        //Finalement, on recommence la boucle while
        while (!indexTrie.empty() and vu_total[indexTrie.top().key]) {
            indexTrie.pop();
        }
        if (!indexTrie.empty()) {
            index = indexTrie.top().key;
            indexTrie.pop();
            vu_total[index] = 1;
        } else {
            break;
        }
    }

    cout << "Fin de la recherche de composantes." << endl;
}