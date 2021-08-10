/**************************************************************************************************
 * Implementation of the TAD Graph
**************************************************************************************************/

#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED
#include "Node.h"
#include <fstream>
#include <stack>
#include <list>

using namespace std;

class Graph{

    //Atributes
private:
    int order;
    int number_edges;
    bool directed;
    bool weighted_edge;
    bool weighted_node;
    Node* first_node;
    Node* last_node;

public:
    //Constructor
    Graph(int order, bool directed, bool weighted_edge, bool weighted_node);
    //Destructor
    ~Graph();
    //Getters
    int getOrder();
    int getNumberEdges();
    bool getDirected();
    bool getWeightedEdge();
    bool getWeightedNode();
    Node* getFirstNode();
    Node* getLastNode();
    //Other methods
    void insertNode(int id);
    void insertEdge(int id, int target_id, float weight);
    void removeNode(int id);
    bool searchNode(int id);
    Node* getNode(int id);

    //methods phase1
    void topologicalSorting();
    void depthFirstSearch(int id, ofstream& output_file);
    Graph* getVertexInduced(int* listIdNodes);
    Graph* agmKruskal(ofstream &output_file);
    Graph* agmPrim(ofstream& output_file);
    float floydMarshall(int idSource, int idTarget, ofstream &output_file);
    float dijkstra(int idSource, int idTarget, ofstream &output_file);

    //methods phase1
    float greed();
    float greedRandom();
    float greedRactiveRandom();

    void printFTD(Node* node, ofstream &output_file, string type);
    Graph* inverseGraph();

private:
    //Auxiliar methods
    void printFTDRec(Node* nodeFTD, int id, ofstream &output_file);
    vector<int> dijkstraDist(int s, vector<int>& path);
    void printPath(vector<int> path, int i, int s, ofstream &output_file, bool &pathFound);
    void printPathFloyd(vector<int>& path, ofstream &output_file);
    int minKey(int key[], bool mstSet[], int V);

};

#endif // GRAPH_H_INCLUDED
