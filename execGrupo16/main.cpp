#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <utility>
#include <tuple>
#include <iomanip>
#include <stdlib.h>
#include <chrono>
#include "Graph.h"
#include "Node.h"

using namespace std;

Graph *leitura(ifstream &input_file, int directed, int weightedEdge, int weightedNode) {

    //Variáveis para auxiliar na criação dos nós no Grafo
    int idNodeSource;
    int idNodeTarget;
    int order;

    //Pegando a ordem do grafo
    input_file >> order;

    //Criando objeto grafo
    Graph *graph = new Graph(order, directed, weightedEdge, weightedNode);

    //Leitura de arquivo

    if (!graph->getWeightedEdge() && !graph->getWeightedNode()) {

        while (input_file >> idNodeSource >> idNodeTarget) {

            graph->insertEdge(idNodeSource, idNodeTarget, 0);

        }

    } else if (graph->getWeightedEdge() && !graph->getWeightedNode()) {

        float edgeWeight;

        while (input_file >> idNodeSource >> idNodeTarget >> edgeWeight) {

            graph->insertEdge(idNodeSource, idNodeTarget, edgeWeight);

        }

    } else if (graph->getWeightedNode() && !graph->getWeightedEdge()) {

        float nodeSourceWeight, nodeTargetWeight;

        while (input_file >> idNodeSource >> nodeSourceWeight >> idNodeTarget >> nodeTargetWeight) {

            graph->insertEdge(idNodeSource, idNodeTarget, 0);
            graph->getNode(idNodeSource)->setWeight(nodeSourceWeight);
            graph->getNode(idNodeTarget)->setWeight(nodeTargetWeight);

        }

    } else if (graph->getWeightedNode() && graph->getWeightedEdge()) {

        float nodeSourceWeight, nodeTargetWeight, edgeWeight;

        while (input_file >> idNodeSource >> nodeSourceWeight >> idNodeTarget >> nodeTargetWeight) {

            graph->insertEdge(idNodeSource, idNodeTarget, edgeWeight);
            graph->getNode(idNodeSource)->setWeight(nodeSourceWeight);
            graph->getNode(idNodeTarget)->setWeight(nodeTargetWeight);

        }

    }

    return graph;
}

//Graph *leituraInstancia(ifstream &input_file, int directed, int weightedEdge, int weightedNode) {
//
//
//    //Variáveis para auxiliar na criação dos nós no Grafo
//    int idNodeSource;
//    int idNodeTarget;
//    int order;
//    int numEdges;
//
//
//    //Pegando a ordem do grafo
//    input_file >> order;
//
//    //Criando objeto grafo
//    Graph *graph = new Graph(order, directed, weightedEdge, weightedNode);
//
//    //Leitura de arquivo
//    while (input_file >> idNodeSource >> idNodeTarget) {
//
//        graph->insertEdge(idNodeSource, idNodeTarget, 0);
//
//    }
//
//    return graph;
//}

void fechoTransitivoDireto(Graph *graph, ofstream &output_file) {
    int id;

    if (!graph->getDirected()) {
        cout << "Graph is not directed" << endl;
    } else {
        cout << "Node id: " << endl;
        cin >> id;
        Node *node = graph->getNode(id);
        if (node == nullptr) {
            cout << "Node not found" << endl;
        } else {
            graph->printFTD(node, output_file, "Direto");
        }
    }

}

void fechoTransitivoIndireto(Graph *graph, ofstream &output_file) {
    int id;

    if (!graph->getDirected()) {
        cout << "Graph is not directed" << endl;
    } else {
        Graph* newGraph = graph->inverseGraph();
        cout << "Node id: " << endl;
        cin >> id;
        Node *node = newGraph->getNode(id);
        if (node == nullptr) {
            cout << "Node not found" << endl;
        } else {
            newGraph->printFTD(node, output_file, "Indireto");
        }
    }
}

void caminhoMinimoDijkstra(Graph *graph, ofstream &output_file) {
    int id1, id2;

    if (!graph->getWeightedEdge()) {
        cout << "Graph is not weighted" << endl;
    } else {
        cout << "Source node id: " << endl;
        cin >> id1;

        cout << "Destination node id: " << endl;
        cin >> id2;

        graph->dijkstra(id1, id2, output_file);
    }
}

void caminhoMinimoFloyd(Graph *graph, ofstream &output_file) {
    int id1, id2;

    if (!graph->getDirected()) {
        cout << "Graph is not directed" << endl;
    } else {
        if (!graph->getWeightedEdge()) {
            cout << "Graph is not weighted" << endl;
        } else {
            cout << "Source node id: " << endl;
            cin >> id1;

            cout << "Destination node id: " << endl;
            cin >> id2;

            graph->floydMarshall(id1, id2, output_file);
        }
    }
}

void arvoreGeradoraMinimaPrim(Graph *graph, ofstream &output_file) {
    if (graph->getDirected()) {
        cout << "Graph is directed" << endl;
    } else {
        if (!graph->getWeightedEdge()) {
            cout << "Graph is not weighted" << endl;
        } else {
            graph->agmPrim(output_file);
        }
    }
}

void arvoreGeradoraMinimaKruskal(Graph *graph, ofstream &output_file) {
    if (!graph->getDirected()) {
        cout << "Graph is directed" << endl;
    } else {
        if (!graph->getWeightedEdge()) {
            cout << "Graph is not weighted" << endl;
        } else {
            graph->agmKruskal(output_file);
        }
    }
}

void buscaProfundidade(Graph *graph, ofstream &output_file) {
    int id1;
    if (!graph->getDirected()) {
        cout << "Graph is directed" << endl;
    } else {
        if (graph->getWeightedEdge()) {
            cout << "Graph is weighted" << endl;
        } else {
            cout << "Source node id: " << endl;
            cin >> id1;
            graph->depthFirstSearch(id1, output_file);
        }
    }
}

void ordenacaoTopologica(Graph *graph) {
    int id1;
    if (!graph->getDirected()) {
        cout << "Graph is directed" << endl;
    } else {
        if (graph->getWeightedEdge()) {
            cout << "Graph is weighted" << endl;
        } else {
            graph->topologicalSorting();
        }
    }
}

int menu() {

    int selecao;

    cout << "MENU" << endl;
    cout << "----" << endl;
    cout << "[1] Fecho transitivo direto de um vértice de um grafo direcionado" << endl;
    cout << "[2] Fecho transitivo indireto de um vértice de um grafo direcionado" << endl;
    cout << "[3] Caminho Mínimo entre dois vértices - Dijkstra" << endl;
    cout << "[4] Caminho Mínimo entre dois vértices - Floyd" << endl;
    cout << "[5] Árvore Geradora Mínima de Prim" << endl;
    cout << "[6] Árvore Geradora Mínima de Kruskal" << endl;
    cout << "[7] Imprimir caminhamento em profundidade" << endl;
    cout << "[8] Imprimir ordenacao topológica" << endl;
    cout << "[0] Sair" << endl;

    cin >> selecao;

    return selecao;

}

void selecionar(int selecao, Graph *graph, ofstream &output_file) {

    switch (selecao) {

        //Fecho transitivo direto de um vértice de um grafo direcionado
        case 1: {
            fechoTransitivoDireto(graph, output_file);
            break;
        }

        //Fecho transitivo indireto de um vértice de um grafo direcionado
        case 2: {
            fechoTransitivoIndireto(graph, output_file);
            break;
        }
        //Caminho mínimo entre dois vértices usando Dijkstra;
        case 3: {
            caminhoMinimoDijkstra(graph, output_file);
            break;
        }

        //Caminho mínimo entre dois vértices usando Floyd;
        case 4: {
            caminhoMinimoFloyd(graph, output_file);
            break;
        }

        //AGM Prim;
        case 5: {
            arvoreGeradoraMinimaPrim(graph, output_file);
            break;
        }

        //AGM - Kruskal;
        case 6: {
            arvoreGeradoraMinimaKruskal(graph, output_file);
            break;
        }

        //Busca em profundidade;
        case 7: {
            buscaProfundidade(graph, output_file);
            break;
        }
        //Ordenação Topologica;
        case 8: {
            ordenacaoTopologica(graph);
            break;
        }

        case 0: {
            break;
        }

        default: {
            cout << "Error!!! invalid option!!" << endl;
        }

    }
}

int mainMenu(ofstream &output_file, Graph *graph) {

    int selecao = 1;

    while (selecao != 0) {
        // system("clear");
        selecao = menu();

        if (output_file.is_open())
            selecionar(selecao, graph, output_file);

        else
            cout << "Unable to open the output_file" << endl;

        output_file << endl;

    }

    return 0;
}


int main(int argc, char const *argv[]) {

    // todos os nós dos grafos inicializam de zero e sao numerados sequencialmente

//    argc = 6;
//    argv[0] = "/Users/juliaaraujo/CLionProjects/Grafo";
//    argv[1] = "/Users/juliaaraujo/CLionProjects/entrada/nao_ponderados/grafo_topologico.txt";
//    argv[2] = "/Users/juliaaraujo/CLionProjects/saida/nao_ponderados/saida_dfs.txt";
//    argv[3] = "1"; // 0 = Nao Direcionado / 1 = Direcionado
//    argv[4] = "0"; // 0 = Nao Ponderado / 1 = Ponderado (arestas)
//    argv[5] = "0"; // 0 = Nao Ponderado / 1 = Ponderado (nos)

    //Verificação se todos os parâmetros do programa foram entrados
    if (argc != 6) {
        cout << "ERROR: Expecting: ./<program_name> <input_file> <output_file> <directed> <weighted_edge> <weighted_node> " << endl;
        return 1;
    }

    string program_name(argv[0]);
    string input_file_name(argv[1]);

    string instance;
    if (input_file_name.find("v") <= input_file_name.size()) {
        string instance = input_file_name.substr(input_file_name.find("v"));
        cout << "Running " << program_name << " with instance " << instance << " ... " << endl;
    }


    //Abrindo arquivo de entrada
    ifstream input_file;
    ofstream output_file;
    input_file.open(argv[1], ios::in);
    output_file.open(argv[2], ios::out | ios::trunc);

    Graph *graph;

    if (input_file.fail()) {
        cout << "Unable to open " << argv[1] << endl;
    }
    else {
        graph = leitura(input_file, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
        mainMenu(output_file, graph);
    }

    //Fechando arquivo de entrada
    input_file.close();

    //Fechando arquivo de saída
    output_file.close();


    return 0;
}


