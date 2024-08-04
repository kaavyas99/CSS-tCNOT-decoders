#include "BasicFunctions.h"


#ifndef DIST_FNC_H
#define DIST_FNC_H

using namespace std;

// A structure to represent a node in adjacency list
struct AdjListNode
{
    int dest;
    int weight;
    double prob;
    AdjListNode* next;
};

// A structure to represent an adjacency list
struct AdjList
{
    AdjListNode *head; // Pointer to head node of list

};

// A structure to represent a graph, an array of adjacency lists.
// Size of array will be V (number of vertices in graph)
struct Graph
{
    int V;
    AdjList* array;

    ~Graph() 
    {
        if (array != NULL) 
        {
            for (int i = 0; i < V; ++i) {
                AdjListNode* current = array[i].head;
                while (current != NULL) {
                    AdjListNode* next = current->next;
                    delete current;
                    current = next;
                }
            }
            delete[] array;
        }
    }


};

AdjListNode* newAdjListNode(int dest, double p);
Graph* createGraph(int V);
void printGraph(Graph* graph);
void printArr(int dist[], int n);
vector<int>  dijkstra(Graph* graph, int src);
int getWeight(double prob);

#endif