#include "Graph.h"
#include "Node.h"
#include "Edge.h"
#include <iostream>
#include <fstream>
#include <stack>
#include <queue>
#include <list>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <float.h>
#include <iomanip>
#include <unordered_set>
#include <map>


#define infi 1000000000
using namespace std;

/**************************************************************************************************
 * Defining the Graph's methods
**************************************************************************************************/

// Constructor
Graph::Graph(int order, bool directed, bool weighted_edge, bool weighted_node)
{

    this->order = order;
    this->directed = directed;
    this->weighted_edge = weighted_edge;
    this->weighted_node = weighted_node;
    this->first_node = this->last_node = nullptr;
    this->number_edges = 0;
}

// Destructor
Graph::~Graph()
{

    Node *next_node = this->first_node;

    while (next_node != nullptr)
    {

        next_node->removeAllEdges();
        Node *aux_node = next_node->getNextNode();
        delete next_node;
        next_node = aux_node;
    }
}

// Getters
int Graph::getOrder()
{

   return this->order;
}
int Graph::getNumberEdges()
{

    return this->number_edges;
}
//Function that verifies if the graph is directed
//Função que verifica se o grafo é direcionado
bool Graph::getDirected()
{

    return this->directed;
}
//Function that verifies if the graph is weighted at the edges
//Função que verifica se o grafo é ponderado nas arestas
bool Graph::getWeightedEdge()
{

    return this->weighted_edge;
}

//Function that verifies if the graph is weighted at the nodes
//Função que verifica se o grafo é ponderado nos nós
bool Graph::getWeightedNode()
{

    return this->weighted_node;
}


Node *Graph::getFirstNode()
{

    return this->first_node;
}

Node *Graph::getLastNode()
{

    return this->last_node;
}

// Other methods
//Outros metodos
/*
    The outdegree attribute of nodes is used as a counter for the number of edges in the graph.
    This allows the correct updating of the numbers of edges in the graph being directed or not.


    O atributo outdegree de nós é usado como um contador para o número de arestas no grafo.
    Isso permite a atualização correta do número de arestas no grafo sendo direcionado ou não.
*/
void Graph::insertNode(int id)
{

    if(this->first_node != nullptr){
        // Allocating the new node and keeping the integrity of the node list
        //Alocando o novo nó mantendo a integridade da lista de nós
        Node* node = new Node(id);
        this->last_node->setNextNode(node);
        this->last_node = node;

    }
    else{
        // Allocating the new node and keeping the integrity of the node list
        //Alocando o novo nó mantendo a integridade da lista de nós
        this->first_node = new Node(id);
        this->last_node = this->first_node;

    }

}

void Graph::insertEdge(int id, int target_id, float weight)
{
    if (!this->searchNode(id)) {
        this->insertNode(id);
    }
    Node* node1 = this->getNode(id);

    if (!this->searchNode(target_id)) {
        this->insertNode(target_id);
    }
    Node* node2 = this->getNode(target_id);

    node1->insertEdge(target_id, weight);
    this->number_edges += 1;

    if (this->directed == 0) {
        node2->insertEdge(id, weight);
        this->number_edges += 1;
    }

}

void Graph::removeNode(int id){

}

bool Graph::searchNode(int id)
{
    Node *next_node = this->first_node;

    while (next_node != nullptr)
    {
        if (next_node->getId() == id) {
            return true;
        }
        next_node = next_node->getNextNode();
    }
    return false;
}

Node *Graph::getNode(int id)
{
    Node *next_node = this->first_node;

    while (next_node != nullptr)
    {
        if (next_node->getId() == id) {
            return next_node;
        }
        next_node = next_node->getNextNode();;
    }
    return nullptr;
}



void Graph::printPathFloyd(vector<int>& path, ofstream &output_file)
{
    int n = path.size();
    for (int i = 0; i < n - 1; i++) {
        cout << path[i] << " -> ";
        output_file << path[i] << " -> ";
    }
    cout << path[n - 1] << endl;
    output_file << path[n - 1] << ";" << endl;
}

float Graph::floydMarshall(int idSource, int idTarget, ofstream &output_file){
    int graphSize = this->order;
    int dis[graphSize][graphSize];
    int next[graphSize][graphSize];

    int matrix[graphSize][graphSize];

    for (int i=0; i < graphSize; i++) {
        for (int j=0; j < graphSize; j++) {
            if (i == j) {
                matrix[i][j] = 0;
            } else {
                matrix[i][j] = infi;
            }
        }
    }

    Node *next_node = this->first_node;
    while (next_node != nullptr)
    {
        Edge* next_edge = next_node->getFirstEdge();
        while(next_edge != nullptr){
            matrix [next_node->getId()][next_edge->getTargetId()] = next_edge->getWeight();
            next_edge = next_edge->getNextEdge();
        }
        next_node = next_node->getNextNode();
    }

    // initialise
    //inicializando

    for (int i = 0; i < graphSize; i++) {
        for (int j = 0; j < graphSize; j++) {
            dis[i][j] = matrix[i][j];

            // No edge between node i and j
            if (matrix[i][j] == infi)
                next[i][j] = -1;
            else
                next[i][j] = j;
        }
    }

    // Standard Floyd Warshall Algorithm
    //Aloritmo de Floyd Warshall

    for (int k = 0; k < graphSize; k++) {
        for (int i = 0; i < graphSize; i++) {
            for (int j = 0; j < graphSize; j++) {

                // We cannot travel through edge that doesn't exist
                // Não podemos percorrer as arestas que não existem
                if (dis[i][k] == infi || dis[k][j] == infi)
                    continue;

                if (dis[i][j] > dis[i][k] + dis[k][j]) {
                    dis[i][j] = dis[i][k] + dis[k][j];
                    next[i][j] = next[i][k];
                }
            }
        }
    }

    // Construct path
    //Construção do caminho

    // If there's no path between node u and v, simply return an empty array
    //Se não houver um caminho entre os nós u e v, simplismente retornar um array vazio
    if (next[idSource][idTarget] == -1)
        return {};

    // Storing the path in a vector
    //Armazenando o caminho em um vetor
    vector<int> path = { idSource };
    while (idSource != idTarget) {
        idSource = next[idSource][idTarget];
        path.push_back(idSource);
    }

    // Print path
    //Imprimir o caminho

    output_file << "digraph caminhoMinimoFloyd {" << endl;
    output_file << "   ";
    this->printPathFloyd(path, output_file);
    output_file << "}" << endl;

}

// Caminho mínimo de Dijkstra

vector<int> Graph::dijkstraDist(int s, vector<int>& path)
{
    // Stores distance of each vertex from source vertex
    //Armazena a distancia de cada vertice do vertice de origem
    vector<int> dist(this->order);

    // Boolean array that shows whether the vertex 'i' is visited or not
    //Array booleano que mostra se o vértice 'i' é visitado ou não
    bool visited[this->order];
    for (int i = 0; i < this->order; i++) {
        visited[i] = false;
        path[i] = -1;
        dist[i] = infi;
    }
    dist[s] = 0;
    path[s] = -1;
    int current = s;

    // Set of vertices that has a parent (one or more) marked as visited
    // Conjunto de vértices que tem um pai (um ou mais) marcado como visitado

    unordered_set<int> sett;
    while (true) {

        // Mark current as visited
        // Marcar atual como visitado

        visited[current] = true;

        Edge* next_edge = this->getNode(current)->getFirstEdge();

        while(next_edge != nullptr){

            int v = next_edge->getTargetId();

            if (visited[v]) {
                next_edge = next_edge->getNextEdge();
                continue;
            }

            // Inserting into the visited vertex
            // Inserindo no vértice visitado

            sett.insert(v);
            int alt = dist[current] + next_edge->getWeight();

            // Condition to check the distance is correct and update it if it is minimum from the previous computed distance
            // Condição para verificar se a distância está correta e atualizá-la se for mínima em relação à distância calculada anteriormente

            if (alt < dist[v]) {
                dist[v] = alt;
                path[v] = current;
            }
            next_edge = next_edge->getNextEdge();
        }

        sett.erase(current);
        if (sett.empty())
            break;

        int minDist = infi;
        int index = 0;

        // Loop to update the distance of the vertices of the graph
        // Loop para atualizar a distância dos vértices do gráfico
        for (int a: sett) {
            if (dist[a] < minDist) {
                minDist = dist[a];
                index = a;
            }
        }
        current = index;
    }
    return dist;
}

// Function to print the path from the source vertex to the destination vertex
// Função para imprimir o caminho do vértice de origem ao vértice de destino

void Graph::printPath(vector<int> path, int i, int s, ofstream &output_file, bool &pathFound)
{
    if (i != s) {

        // Condition to check if there is no path between the vertices
        // Condição para verificar se não há caminho entre os vértices

        if (path[i] == -1) {
            cout << "Path not found!!" << endl;
            pathFound = false;
            return;
        }
        printPath(path, path[i], s, output_file, pathFound);
        cout << path[i] << " -> ";
        output_file << path[i] << " -> ";
    }
}


float Graph::dijkstra(int idSource, int idTarget, ofstream &output_file){

    if (idSource < 0 || idSource > this->order - 1 || idTarget < 0 || idTarget > this->order - 1) {
        cout << "Path not found!!" << endl;
        return 0;
    }

    output_file << "digraph caminhoMinimoDijkstra {" << endl;
    output_file << "   ";

    vector<int> path(this->order);
    vector<int> dist = dijkstraDist(idSource, path);

    bool pathFound = true;
    printPath(path, idTarget, idSource, output_file, pathFound);

    if (pathFound) {
        cout << idTarget << endl;
        output_file << idTarget << ";" << endl;
    }
    output_file << "}" << endl;

}

// Ordenacao Topologica

// Class to represent a graph
//Classe para representar um grafo
class TopologicalGraph {
    // No. of vertices'
    //Numero de vertices
    int V;

    // Pointer to an array containing adjacency listsList
    // Ponteiro para um array contendo listas de adjacência

    list<int>* adj;

    // A function used by topologicalSort
    // Uma função usada por topologicalSort

    void topologicalSortUtil(int v, bool visited[],
                             stack<int>& Stack);

public:
    // Constructor
    //Construtor
    TopologicalGraph(int V);

    // function to add an edge to graph
    //Função para acrescentar uma aresta ao grafo
    void addEdge(int v, int w);

    // prints a Topological Sort of the complete graph
    // imprime uma classificação topológica do gráfico completo

    void topologicalSort();
};

TopologicalGraph::TopologicalGraph(int V)
{
    this->V = V;
    adj = new list<int>[V];
}

void TopologicalGraph::addEdge(int v, int w)
{
    // Add w to v’s list.
    // Adicione w à lista de v.

    adj[v].push_back(w);
}

// A recursive function used by topologicalSort
// Uma função recursiva usada por topologicalSort

void TopologicalGraph::topologicalSortUtil(int v, bool visited[],
                                stack<int>& Stack)
{
    // Mark the current node as visited.
    // Marque o nó atual como visitado.

    visited[v] = true;

    // Recur for all the vertices adjacent to this vertex
    // Recorre para todos os vértices adjacentes a este vértice

    list<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i])
            topologicalSortUtil(*i, visited, Stack);

    // Push current vertex to stack which stores result
    // Empurra o vértice atual para empilhar o que armazena o resultado

    Stack.push(v);
}

// The function to do Topological Sort.
// A função para fazer a classificação topológica.

// It uses recursive topologicalSortUtil()
// Ele usa topologicalSortUtil () recursivo

void TopologicalGraph::topologicalSort()
{
    stack<int> Stack;

    // Mark all the vertices as not visited
    // Marque todos os vértices como não visitados

    bool* visited = new bool[V];
    for (int i = 0; i < V; i++)
        visited[i] = false;

    // Call the recursive helper function to store Topological
    // Sort starting from all vertices one by one

    // Chame a função auxiliar recursiva para armazenar topológico
    // Classificar começando por todos vértices um por um

    for (int i = 0; i < V; i++)
        if (visited[i] == false)
            topologicalSortUtil(i, visited, Stack);

    // Print contents of stack
    // Imprime o conteúdo da pilha

    while (Stack.empty() == false) {
        cout << Stack.top() << " ";
        Stack.pop();
    }
}


//function that prints a topological sorting
// função que imprime uma classificação topológica

void Graph::topologicalSorting(){

    TopologicalGraph g(this->getOrder());

    Node *next_node = this->first_node;
    while (next_node != nullptr)
    {
        Edge* next_edge = next_node->getFirstEdge();
        while(next_edge != nullptr){

            g.addEdge(next_node->getId(), next_edge->getTargetId());

            next_edge = next_edge->getNextEdge();
        }
        next_node = next_node->getNextNode();
    }

    g.topologicalSort();

    cout << endl;
}


// DTS

class DTSGraph
{
public:
    map<int, bool> visited;
    map<int, list<int>> adj;

    // function to add an edge to graph
    // função para adicionar uma borda ao gráfico

    void addEdge(int v, int w);

    // DFS traversal of the vertices reachable from v
    // Traversal de DFS dos vértices alcançáveis ​​de v

    void DFS(int v, ofstream& output_file);
};

void DTSGraph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    // Add w to v’s list.
    // Adicione w à lista de v.

}

void DTSGraph::DFS(int v, ofstream& output_file)
{
    // Mark the current node as visited and print it
    // Marque o nó atual como visitado e imprima
    visited[v] = true;
    cout << v ;
    output_file << v;

    // Recur for all the vertices adjacent to this vertex
    // Recorre para todos os vértices adjacentes a este vértice

    list<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i]) {
            cout << " -> ";
            output_file << " -> ";
            DFS(*i, output_file);
        }
}


void Graph::depthFirstSearch(int id, ofstream& output_file){
    DTSGraph g;

    Node *next_node = this->first_node;
    while (next_node != nullptr)
    {
        Edge* next_edge = next_node->getFirstEdge();
        while(next_edge != nullptr){

            g.addEdge(next_node->getId(), next_edge->getTargetId());

            next_edge = next_edge->getNextEdge();
        }
        next_node = next_node->getNextNode();
    }

    output_file << "digraph caminhoProfundidade {" << endl;
    output_file << "   ";

    g.DFS(id, output_file);

    cout << endl;

    output_file << endl;
    output_file << "}" << endl;

}


// Arvore Geradora Minima de Kruskal

// a structure to represent a weighted edge in graph
// uma estrutura para representar um borda ponderada no gráfico
class EdgeKruskal {
public:
    int src, dest, weight;
};

// a structure to represent a connected, undirected and weighted graph
// uma estrutura para representar um conectado, gráfico não direcionado e ponderado
class GraphKruskal {
public:

    // V-> Number of vertices, E-> Number of edges
    // V-> Número de vértices, E-> Número de arestas
    int V, E;

    // graph is represented as an array of edges.
    // Since the graph is undirected, the edge from src to dest is also edge from dest to src. Both are counted as 1 edge here.
    // gráfico é representado como uma matriz de arestas.

    // Uma vez que o gráfico não é direcionado, a aresta de src para dest também é aresta de dest para src. Ambos são contados como 1 aresta aqui.
    EdgeKruskal* edge;
};

// Creates a graph with V vertices and E edges
// Cria um gráfico com vértices V e arestas E

GraphKruskal* createGraph(int V, int E)
{
    GraphKruskal* graph = new GraphKruskal;
    graph->V = V;
    graph->E = E;

    graph->edge = new EdgeKruskal[E];

    return graph;
}

// A structure to represent a subset for union-find
// Uma estrutura para representar um subconjunto para union-find

class subset {
public:
    int parent;
    int rank;
};

// A utility function to find set of an element i
// (uses path compression technique)


// Uma função de utilidade para encontrar o conjunto de um elemento i
// (usa técnica de compressão de caminho)
int find(subset subsets[], int i)
{
    // find root and make root as parent of i
    // (path compression)

    // encontre a raiz e faça a raiz como pai de i
    // (compressão de caminho)
    if (subsets[i].parent != i)
        subsets[i].parent = find(subsets, subsets[i].parent);

    return subsets[i].parent;
}

// A function that does union of two sets of x and y (uses union by rank)

// Uma função que faz a união de dois conjuntos de x e y (usa união por classificação)
void Union(subset subsets[], int x, int y)
{
    int xroot = find(subsets, x);
    int yroot = find(subsets, y);

    // Attach smaller rank tree under root of high rank tree (Union by Rank)

    // Anexe a árvore de classificação menor sob a raiz de alta árvore de classificação (União por Classificação)
    if (subsets[xroot].rank < subsets[yroot].rank)
        subsets[xroot].parent = yroot;
    else if (subsets[xroot].rank > subsets[yroot].rank)
        subsets[yroot].parent = xroot;

        // If ranks are same, then make one as root and increment its rank by one

        // Se as classificações forem iguais, faça uma como root e incrementa sua classificação em um
    else {
        subsets[yroot].parent = xroot;
        subsets[xroot].rank++;
    }
}

// Compare two edges according to their weights.
// Used in qsort() for sorting an array of edges

// Compare duas arestas de acordo com seus pesos.
// Usado em qsort () para classificar uma matriz de arestas
int myComp(const void* a, const void* b)
{
    EdgeKruskal* a1 = (EdgeKruskal*)a;
    EdgeKruskal* b1 = (EdgeKruskal*)b;
    return a1->weight > b1->weight;
}

// The main function to construct MST using Kruskal's algorithm
// A principal função para construir MST usando Kruskal's algoritmo
void KruskalMST(GraphKruskal* graph, ofstream &output_file)
{
    int V = graph->V;
    EdgeKruskal result[V]; // Tnis will store the resultant MST // Tnis armazenará o MST resultante
    int e = 0; // An index variable, used for result[] // Uma variável de índice, usada para o resultado []
    int i = 0; // An index variable, used for sorted edges // Uma variável de índice, usada para bordas classificadas

    // Step 1: Sort all the edges in non-decreasing
    // order of their weight. If we are not allowed to
    // change the given graph, we can create a copy of
    // array of edges

    // Etapa 1: classificar todas as arestas em não decrescentes
    // ordem de seu peso. Se não temos permissão para
    // alterar o gráfico fornecido, podemos criar uma cópia de
    // array de bordas
    qsort(graph->edge, graph->E, sizeof(graph->edge[0]),
          myComp);

    // Allocate memory for creating V ssubsets
    // Alocar memória para criar subsets V
    subset* subsets = new subset[(V * sizeof(subset))];

    // Create V subsets with single elements
    // Cria subconjuntos V com elementos únicos

    for (int v = 0; v < V; ++v)
    {
        subsets[v].parent = v;
        subsets[v].rank = 0;
    }

    // Number of edges to be taken is equal to V-1
    // O número de arestas a serem tomadas é igual a V-1

    while (e < V - 1 && i < graph->E)
    {
        // Step 2: Pick the smallest edge. And increment
        // the index for next iteration
        // Etapa 2: escolha a menor aresta. E incremento
        // o índice para a próxima iteração
        EdgeKruskal next_edge = graph->edge[i++];

        int x = find(subsets, next_edge.src);
        int y = find(subsets, next_edge.dest);

        // If including this edge does't cause cycle,
        // include it in result and increment the index
        // of result for next edge

        // Se incluir esta aresta não causa ciclo,
        // inclua no resultado e incremente o índice
        // do resultado para a próxima borda
        if (x != y) {
            result[e++] = next_edge;
            Union(subsets, x, y);
        }
        // Else discard the next_edge
        // Caso contrário, descarte o next_edge

    }

    // print the contents of result[] to display the
    // built MST
    // imprime o conteúdo do resultado [] para exibir o
    // construído MST


    cout<<"Edge \tWeight\n";
    output_file << "graph arvoreGeradoraMinimaKruskal {" << endl;

    for (i = 0; i < e; ++i)
    {
        cout << result[i].src << " - " << result[i].dest << " \t" << result[i].weight << endl;
        output_file << "   " << result[i].src << " -- " << result[i].dest << ";" << endl;
    }

    output_file << "}" << endl;
}

Graph* Graph::agmKruskal(ofstream &output_file){

    int nodeNumber = 0;

    GraphKruskal* graph = createGraph(this->order, this->number_edges);

    Node *next_node = this->first_node;
    while (next_node != nullptr)
    {
        Edge* next_edge = next_node->getFirstEdge();
        while(next_edge != nullptr){

            graph->edge[nodeNumber].src = next_node->getId();
            graph->edge[nodeNumber].dest = next_edge->getTargetId();
            graph->edge[nodeNumber].weight = next_edge->getWeight();
            nodeNumber++;

            next_edge = next_edge->getNextEdge();
        }
        next_node = next_node->getNextNode();
    }

    // Function call
    //Chamada de função
    KruskalMST(graph, output_file);
}


// Arvore Geradora Minima Prim

// A utility function to find the vertex with
// minimum key value, from the set of vertices
// not yet included in MST


// Uma função de utilidade para encontrar o vértice com
// valor-chave mínimo, do conjunto de vértices
// ainda não incluído no MST

int Graph::minKey(int key[], bool mstSet[], int V)
{
    // Initialize min value
    int min = INT_MAX, min_index;

    for (int v = 0; v < V; v++)
        if (mstSet[v] == false && key[v] < min)
            min = key[v], min_index = v;

    return min_index;
}


Graph* Graph::agmPrim(ofstream& output_file){
    int graphSize = this->order;
    int matrix[graphSize][graphSize];

    for (int i=0; i < graphSize; i++) {
        for (int j=0; j < graphSize; j++) {
            if (i == j) {
                matrix[i][j] = 0;
            } else {
                matrix[i][j] = infi;
            }
        }
    }

    Node *next_node = this->first_node;
    while (next_node != nullptr)
    {
        Edge* next_edge = next_node->getFirstEdge();
        while(next_edge != nullptr){
            matrix [next_node->getId()][next_edge->getTargetId()] = next_edge->getWeight();
            next_edge = next_edge->getNextEdge();
        }
        next_node = next_node->getNextNode();
    }


    // Prim MST

    // Array to store constructed MST
    // Matriz para armazenar MST construído

    int parent[graphSize];

    // Key values used to pick minimum weight edge in cut
    // Valores-chave usados ​​para escolher a borda de peso mínimo no corte

    int key[graphSize];

    // To represent set of vertices included in MST
    // Para representar o conjunto de vértices incluídos no MST

    bool mstSet[graphSize];

    // Initialize all keys as INFINITE
    // Inicializa todas as chaves como INFINITE

    for (int i = 0; i < graphSize; i++)
        key[i] = INT_MAX, mstSet[i] = false;

    // Always include first 1st vertex in MST.
    // Make key 0 so that this vertex is picked as first vertex.
    // Sempre inclua o primeiro primeiro vértice no MST.
    // Faça a chave 0 para que este vértice seja escolhido como o primeiro vértice.
    key[0] = 0;
    parent[0] = -1; // First node is always root of MST // O primeiro nó é sempre a raiz do MST

    // The MST will have V vertices
    // O MST terá vértices V
    for (int count = 0; count < graphSize - 1; count++)
    {
        // Pick the minimum key vertex from the
        // set of vertices not yet included in MST
        // Escolha o vértice chave mínimo do
        // conjunto de vértices ainda não incluídos no MST
        int u = minKey(key, mstSet, graphSize);

        // Add the picked vertex to the MST Set
        // Adicione o vértice escolhido ao conjunto MST

        mstSet[u] = true;

        // Update key value and parent index of
        // the adjacent vertices of the picked vertex.
        // Consider only those vertices which are not
        // yet included in MST

        // Atualize o valor da chave e o índice pai de
        // os vértices adjacentes do vértice escolhido.
        // Considere apenas os vértices que não são
        // ainda incluído no MST
        for (int v = 0; v < graphSize; v++)

            // graph[u][v] is non zero only for adjacent vertices of m
            // mstSet[v] is false for vertices not yet included in MST
            // Update the key only if graph[u][v] is smaller than key[v]


            // gráfico [u] [v] é diferente de zero apenas para vértices adjacentes de m
            // mstSet [v] é falso para vértices ainda não incluídos no MST
            // Atualize a chave apenas se o gráfico [u] [v] for menor que a chave [v]
            if (matrix[u][v] && mstSet[v] == false && matrix[u][v] < key[v])
                parent[v] = u, key[v] = matrix[u][v];
    }

    // print the constructed MST

    // A utility function to print the
    // constructed MST stored in parent[]
    // imprime o MST construído

    // Uma função de utilidade para imprimir o
    // construído MST armazenado em pai []


    cout<<"Edge \tWeight\n";
    output_file << "graph arvoreGeradoraMinimaPrim {" << endl;

    for (int i = 1; i < graphSize; i++) {
        cout << parent[i] << " - " << i << " \t" << matrix[i][parent[i]] << endl;
        output_file << "   " << parent[i] << " -- " << i << ";" << endl;
    }

    output_file << "}" << endl;
}

void Graph::printFTDRec(Node* nodeFTD, int id, ofstream &output_file) {
    Node *node = this->getNode(id);
    for (Edge *auxEdge = node->getFirstEdge(); auxEdge != nullptr; auxEdge = auxEdge->getNextEdge()) {
        if (!nodeFTD->searchEdge(auxEdge->getTargetId())) {
            nodeFTD->insertEdge(auxEdge->getTargetId(), 0);
            cout << auxEdge->getTargetId() << " ";
            output_file << "   " << auxEdge->getTargetId() << ";" << endl;
            printFTDRec(nodeFTD, auxEdge->getTargetId(), output_file);
        }
    }
}

void Graph::printFTD(Node* node, ofstream &output_file, string type) {
    Node* nodeFTD = new Node(0);

    output_file << "graph fechoTransitivo" << type << " {" << endl;

    printFTDRec(nodeFTD, node->getId(), output_file);

    cout << endl;
    output_file << "}" << endl;
}

Graph* Graph::inverseGraph()
{

    Graph *newGraph = new Graph(this->order, this->directed, this->weighted_edge, this->weighted_node);

    Node *next_node = this->first_node;

    while (next_node != nullptr) {
        for (Edge *aux = next_node->getFirstEdge(); aux != nullptr; aux = aux->getNextEdge()) {
            newGraph->insertEdge(aux->getTargetId(), next_node->getId(), 0);
        }

        next_node = next_node->getNextNode();
    }

    return newGraph;
}

