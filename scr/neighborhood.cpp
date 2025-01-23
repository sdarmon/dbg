/*
 *
 */

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <set>
#include "graph.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wstring-compare"
#define MAX 1024


//Fonction permettant d'enregistrer une composante du graphe `G` dans un fichier
void save_comp(Graph &G, vector<int> &vu, string outputPrefix){
    ofstream output;
    output.open(outputPrefix+".nodes");
    for (int i = 0; i < G.N ; i++){
        if (vu[i]){
            output << i << "\t" << G.Vertices[i].label << "\t" << ((vu[i] < 0) ? 0 : vu[i]) << "\t" << ((vu[i] > 1) ? 0 : G.Vertices[i].weight) << endl;
        }
    }
    return;
}

int main(int argc, char** argv) {
    int threshold,dis;
    string prefixOutput;
    if (argc == 4) {
        prefixOutput = "./";
        threshold = 0;
        dis = 1;
    }else if (argc == 6 and (string)argv[4] == "-o"){
        prefixOutput = argv[5];
        threshold = 0;
        dis = 1;
    }else if (argc == 6 and (string)argv[4] == "-c"){
        prefixOutput = "./";
        threshold = atoi(argv[5]);
        dis = 1;
    }else if (argc == 6 and (string)argv[4] == "-d"){
        prefixOutput = "./";
        threshold = 0;
        dis = atoi(argv[5]);
    }else if (argc == 8 and (string)argv[4] == "-o" and (string)argv[6] == "-d"){
        prefixOutput = argv[5];
        dis = atoi(argv[7]);
        threshold = 0;
    }else if (argc == 8 and (string)argv[4] == "-c" and (string)argv[6] == "-d"){
        threshold = atoi(argv[5]);
        dis = atoi(argv[7]);
        prefixOutput = "./";
    }else if (argc == 10 and (string)argv[4] == "-o" and (string)argv[6] == "-c" and (string)argv[8] == "-d"){
        prefixOutput = argv[5];
        threshold = atoi(argv[7]);
        dis = atoi(argv[9]);
    } else {
        cout << "Expected use of this program: \n\n\t" << argv[0]
             << " file.nodes file.edges node_id/input.txt [-o prefixOutput] [-c value] [-d dis] \n"
             << "The option c allows to give a threshold to the weight of visited vertices (default 0) and the option d gives the "
                "depth of neighbor gathered (default 1)" << endl;
        cout << argc << " arguments have been used " ;
        for(int i =0; i<argc ; i++){
            cout << "arg" << i << ": " << argv[i] << "\t" ;
        }
        cout << endl;
        return 0;
    }

    //On charge le graphe
    vector <Edge> E;
    vector <Node> V;

    char *nodesPath = argv[1];
    char *edgesPath = argv[2];
    vector<int> nodes_id;
    vector<int> pronf;
    string arg = argv[3];
    if (arg.find('.') != string::npos){
        char *idPath =argv[3];
        ifstream ids(idPath, std::ios::binary);
        string line;
        while (getline(ids, line)){
            nodes_id.push_back(stoi(line));
            pronf.push_back(1);
        }
        ids.close();
    } else {
        nodes_id.push_back(stoi(argv[3]));
    }
    ifstream edges(edgesPath, std::ios::binary);
    ifstream nodes(nodesPath, std::ios::binary);

    read_edge_file(edges, E);
    read_node_file_weighted(nodes, V);

    char aretFF[3] = {'F', 'F'};
    char aretFR[3] = {'F', 'R'};
    char aretRF[3] = {'R', 'F'};
    char aretRR[3] = {'R', 'R'};
    Graph G(V, E);
    vector<Edge> E2;
    int node ;
    int p ;
    vector<int> vu(G.N,0);
    //On fait un BFS
    while (not nodes_id.empty() ) { //Cas de terminaison, on a terminé le BFS
        node = nodes_id.front();
        nodes_id.erase(nodes_id.begin());
        p = pronf.front();
        pronf.erase(pronf.begin());
        if (vu[node] != 0) { //Cas où le sommet a été vu par le BFS
            continue;
        }
        vu[node] = p;
        //Case where node is in comp but without any neighbors
        if (p == 1 && G.Neighbors(node)->empty()){
            E2.push_back(Edge(node,node,0,aretFF));
        }
        for (vector<Neighbor>::iterator it = G.Neighbors(node)->begin();
             it != G.Neighbors(node)->end(); ++it) {
            //On boucle sur ses voisins
            if (vu[it->val] == 0 and p + 1 <= dis) {
                //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
                nodes_id.push_back(it->val);
                pronf.push_back(p + 1);
                E2.push_back(Edge(node,it->val,0,it->label));
                if (it->label[0] == 'F' and it->label[1] == 'F' ){
                    E2.push_back(Edge(it->val,node,0,aretRR));
                } else if (it->label[0] == 'F' and it->label[1] == 'R' ){
                    E2.push_back(Edge(it->val,node,0,aretFR));
                } else if (it->label[0] == 'R' and it->label[1] == 'R' ){
                    E2.push_back(Edge(it->val,node,0,aretFF));
                } else if (it->label[0] == 'R' and it->label[1] == 'F' ){
                    E2.push_back(Edge(it->val,node,0,aretRF));
                }
            }
            if (vu[it->val] != 0) {
                //Cas où le voisin a déjà été vu
                E2.push_back(Edge(node,it->val,0,it->label));
                if (it->label[0] == 'F' and it->label[1] == 'F' ){
                    E2.push_back(Edge(it->val,node,0,aretRR));
                } else if (it->label[0] == 'F' and it->label[1] == 'R' ){
                    E2.push_back(Edge(it->val,node,0,aretFR));
                } else if (it->label[0] == 'R' and it->label[1] == 'R' ){
                    E2.push_back(Edge(it->val,node,0,aretFF));
                } else if (it->label[0] == 'R' and it->label[1] == 'F' ){
                    E2.push_back(Edge(it->val,node,0,aretRF));
                }
            }
        }
    }

    save_comp(G, vu, prefixOutput);
    ofstream edges2;
    edges2.open(prefixOutput+".edges");
    printEdges(E2,edges2);
}