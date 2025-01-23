/*
 * Ce programme permet de faire une agglomeration des composantes
 * des graphes.
 */

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>
#include <queue>
#include "graph.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wstring-compare"
#define MAX 1024

//Initialisation d'un vecteur avec la valeur -1. Il me semble que vector<int> vec(-1,n) doit fonctionner, à corriger!!
void initVec(vector<int> &vec, int n){
    vec.clear();
    for(int i = 0; i<n; i++){
        vec.push_back(-1);
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
             << " file.nodes file.edges -k kmer compoPrefix nbComp outputPrefix \n" << endl;
        return 0;
    }

    //On charge le graphe
    vector <Edge> E;
    vector <Node> V;

    string nodesPath = argv[1];
    string edgesPath = argv[2];
    string outputPrefix;
    string compoPrefix;
    int nbComp;
    ifstream edges(edgesPath, std::ios::binary);
    ifstream nodes(nodesPath, std::ios::binary);

    read_edge_file(edges, E);
    read_node_file_weighted(nodes, V);

    Graph G(V, E);
    G.kmer = stoi(argv[4]);
    compoPrefix = argv[5];
    nbComp = stoi(argv[6]);
    outputPrefix = argv[7];


    cout << "Graphe chargé et construit" << endl;

    cout << "Début construction graphe aggloméré" << endl;
    vector <Edge> E3;
    vector <Node> V3;
    E3.clear();
    V3.clear();
    vector<int> seen(G.N, 0);
    vector<int> correspondingVertex(G.N, 0);
    //seen est le tableau de correspondance entre les anciens sommets et les nouveaux.


    //Après, on fait un BFS à partir de chaque sommet dans la composante afin d'obtenir les voisins du périmètre.
    set<int> vu2;
    set<int> sonSet;
    typedef map<pair<int,int>,string> dicChem;
    dicChem areteFF;
    dicChem areteFR;
    dicChem areteRF;
    dicChem areteRR;
    int limiteAretes;
    vector<int> sons;
    vector<int> depthSons;
    queue<int> depth;
    vector<string> labelSons;
    queue<string> labels;
    vector<Neighbor*> aretes;
    queue < Neighbor * > aVoir;
    int index = 0;
    int compteurDeBoucle = 0;
    char aretFF[3] = {'F', 'F'};
    char aretFR[3] = {'F', 'R'};
    char aretRF[3] = {'R', 'F'};
    char aretRR[3] = {'R', 'R'};
    vector <int> comp;
    string line;
    string file_name;
    int compt;
    //On boucle sur les composantes
    cout << "Début des calculs des arêtes des comp " << nbComp<<endl;
    for (int component = 0; component<nbComp; component++){
        comp.clear();
        file_name = compoPrefix+to_string(component)+".txt";
        ifstream file(file_name, std::ios::binary);
        while (getline(file, line)){
            compt=stoi(line.substr(0,line.find('\t')));
            comp.push_back(compt);
            seen[compt] = component + 1;
        }
        areteFF.clear();
        areteFR.clear();
        areteRR.clear();
        areteRF.clear();
        sonSet.clear();
        //On boucle sur les sommets de la composante
        cout << "Pré-calcul pour la composante " << compteurDeBoucle << " de taille " << comp.size() << endl;
        compteurDeBoucle++;
        for (vector<int>::iterator it = comp.begin(); it != comp.end(); ++it) {
            if (G.Vertices[(*it)].label.size() < G.kmer) {
                continue;
            }
            //Cas sommet en périphérie
            sons.clear();
            aretes.clear();
            depthSons.clear();
            labelSons.clear();
            while (!aVoir.empty()) { // N'est jamais censé arriver !
                aVoir.pop();
            }
            while (!depth.empty()) { // N'est jamais censé arriver !
                depth.pop();
            }
            while (!labels.empty()) { // N'est jamais censé arriver !
                labels.pop();
            }
            vu2.clear();
            vu2.insert((*it)); //On voit bien le sommet duquel on part

            //On fait un BFS à partir de chaque sommet afin de savoir quels sommets du périmètre sont
            //atteignables à partir de chaque sommet de la composante
            //Cas en partant du forward
            for (vector<Neighbor>::iterator voisin = G.Neighbors((*it))->begin();
                 voisin != G.Neighbors((*it))->end(); ++voisin) {
                if (voisin->label[0] == 'F') {
                    aVoir.push(&(*voisin));
                    depth.push(1);
                    string aux = G.Vertices[(*it)].label;
                    labels.push(aux);
                }
            }
            G.BFS_comp(seen, vu2, aVoir, depth, labels, sons, depthSons, aretes, labelSons);
            for (int i = 0; i < sons.size(); i++) {
                sonSet.insert(sons[i]);
            }

            limiteAretes = aretes.size(); //On fait ça pour garder en mémoire le fait que ce ne sont
            //pas les mêmes arêtes

            //Cas en partant du reverse
            for (vector<Neighbor>::iterator voisin = G.Neighbors((*it))->begin();
                 voisin != G.Neighbors((*it))->end(); ++voisin) {
                if (voisin->label[0] == 'R') {
                    aVoir.push(&(*voisin));
                    depth.push(1);
                    labels.push("");
                }
            }
            G.BFS_comp(seen, vu2, aVoir, depth, labels, sons, depthSons, aretes, labelSons);

            for (int i = limiteAretes; i < sons.size(); i++) {
                sonSet.insert(sons[i]);
            }
            for (int i = 0; i < limiteAretes; i++) {
                for (int j = limiteAretes; j < sons.size(); j++) {
                    string aux = reverse_complement(labelSons[j]) + labelSons[i];
                    if (aretes[i]->label[1] == 'F' and aretes[j]->label[1] == 'F' and
                        not G.neigh(sons[j],sons[i],'R','F')) { //attention, ici la seconde lettre doit être comprise
                        //comme le complément ! (Voir cahier, ce n'est pas forcément évident)
                        if (areteRF.find(make_pair(sons[j], sons[i])) == areteRF.end() or
                            areteRF[make_pair(sons[j], sons[i])].size() > aux.size()) {
                            areteRF[make_pair(sons[j], sons[i])] = aux;
                        }
                    } else if (aretes[i]->label[1] == 'F' and aretes[j]->label[1] == 'R' and
                        not G.neigh(sons[j],sons[i],'F','F')) {
                        if (areteFF.find(make_pair(sons[j], sons[i])) == areteFF.end() or
                            areteFF[make_pair(sons[j], sons[i])].size() > aux.size()) {
                            areteFF[make_pair(sons[j], sons[i])] = aux;
                        }
                    } else if (aretes[i]->label[1] == 'R' and aretes[j]->label[1] == 'R' and
                        not G.neigh(sons[j],sons[i],'F','R')) {
                        if (areteFR.find(make_pair(sons[j], sons[i])) == areteFR.end() or
                            areteFR[make_pair(sons[j], sons[i])].size() > aux.size()) {
                            areteFR[make_pair(sons[j], sons[i])] = aux;
                        }
                    } else if (not G.neigh(sons[j],sons[i],'R','R')){
                        if (areteRR.find(make_pair(sons[j], sons[i])) == areteRR.end() or
                            areteRR[make_pair(sons[j], sons[i])].size() > aux.size()) {
                            areteRR[make_pair(sons[j], sons[i])] = aux;
                        }
                    }
                }
            }

        } //On termine de traiter tous les sommets de la composante
        cout << "BFS sur tous les sommets terminés" << endl;


        //On peut donc passer la construction du graphe. Commençons par les sommets en péri.
        for (set<int>::iterator setIt = sonSet.begin(); setIt != sonSet.end(); ++setIt) {
            V3.push_back(Node(index, G.Vertices[(*setIt)].weight, G.Vertices[(*setIt)].label));
            correspondingVertex[(*setIt)] = index;
            index++;
            seen[(*setIt)] = -compteurDeBoucle;
        }

        //Puis les arêtes au sein de la composante :
        for (dicChem::iterator itDic = areteFF.begin(); itDic != areteFF.end(); ++itDic) {
            V3.push_back(Node(index, 0, itDic->second));
            E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretFF));
            E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFF));
            E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretRR));
            E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretRR));
            index++;
        }
        for (dicChem::iterator itDic = areteFR.begin(); itDic != areteFR.end(); ++itDic) {
            V3.push_back(Node(index, 0, itDic->second));
            E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretFF));
            E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFR));
            E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretFR));
            E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretRR));
            index++;
        }
        for (dicChem::iterator itDic = areteRF.begin(); itDic != areteRF.end(); ++itDic) {
            if (areteFR.find(make_pair(itDic->first.second, itDic->first.first)) == areteFR.end()) {
                //On évite de créer une arête labellisée supplémentaire si son complémentaire est déjà présent
                //Il doit y avoir des optimisations possibles à faire ici
                V3.push_back(Node(index, 0, itDic->second));
                E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretRF));
                E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFF));
                E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretRR));
                E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretRF));
                index++;
            }
        }
        for (dicChem::iterator itDic = areteRR.begin(); itDic != areteRR.end(); ++itDic) {
            if (areteFF.find(make_pair(itDic->first.second, itDic->first.first)) == areteFF.end()){
                //Idem
                V3.push_back(Node(index, 0, itDic->second));
                E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretRF));
                E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFR));
                E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretFR));
                E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretRF));
                index++;
            }
        }
    } //On vient de terminer cette composante !

    //Maintenant que les composantes ont bien été ajouté, on s'occupe des sommets restants
    for (int i = 0; i < G.N; i++) {
        if (seen[i] == 0 and G.Vertices[i].label.size()>=G.kmer) {
            V3.push_back(Node(index, G.Vertices[i].weight, G.Vertices[i].label));
            correspondingVertex[i] = index;
            index++;
        } else if (seen[i] < 0) {
            continue; //Cas déjà traité précédemment
        } else {
            correspondingVertex[i] = -1; //Cas sommet de comp ou de label trop court
        }
    }
    cout << "Sommets restants ajoutés" << endl;
    limiteAretes = E3.size();
    //Maintenant on s'occupe des arêtes
    for (int i = 0; i<G.N ; i++) {
        if (seen[i] == 0 and correspondingVertex[i]>=0){//Cas sommet non vu et valide
            for(vector<Neighbor>::iterator node = G.Neighbors(i)->begin(); node != G.Neighbors(i)->end(); ++node){
                if (correspondingVertex[node->val]>=0){ //Cas voisin valide
                    E3.push_back(Edge(correspondingVertex[i],correspondingVertex[node->val],0,node->label));
                }
            }
        } else if (seen[i] < 0 ){ //Cas sommet de comp en péri
            for (vector<Neighbor>::iterator node = G.Neighbors(i)->begin(); node != G.Neighbors(i)->end(); ++node) {
                if (seen[node->val] <= 0 and correspondingVertex[node->val]>=0) { //Voisin hors comp et valide
                    E3.push_back(Edge(correspondingVertex[i], correspondingVertex[node->val], 0, node->label));
                }
            }
        }
    }
    cout << "Arêtes restantes ajoutées" << endl;

    //Enfin, on enregistre le graphe créé.
    ofstream outputNodes2;
    outputNodes2.open(outputPrefix + "/clean.nodes");
    printVerticesBcalm(V3, outputNodes2);
    outputNodes2.close();

    ofstream outputEdges2;
    outputEdges2.open(outputPrefix + "/clean.edges");
    printEdgesBcalm(E3, outputEdges2);
    outputEdges2.close();

    //On récupère aussi l'abondance qui est nécessaire pour kissplice
    ifstream ab(edgesPath.substr(0,edgesPath.size()-12)+".abundance", std::ios::binary);

    vector <double> A;
    read_abundance_file(ab,A);
    vector<double> A3(V3.size(),0.0);
    for (int i = 0; i<A.size();i++){
        if(correspondingVertex[i] >= 0){
            A3[correspondingVertex[i]]+= A[i];
        }
    }

    ofstream outputAb;
    outputAb.open(outputPrefix + "/clean.abundance");
    printAbundance(A3,outputAb);
    cout << "Graphe enregistré" << endl;

}