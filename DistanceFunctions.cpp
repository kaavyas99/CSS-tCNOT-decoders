#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <limits.h>
#include <vector>
#include <set>
#include "BasicFunctions.h"
#include "DistanceFunctions.hpp"
using namespace std;

/***************************************************************************************
*    C++ program for Dijkstra's  shortest path algorithm for adjacency list representation of graph
*    Uses open source code available online
***************************************************************************************/



// A utility function to create a new adjacency list node
AdjListNode* newAdjListNode(int dest, double prob)
{
    AdjListNode* newNode = new AdjListNode();
    newNode->dest = dest;
    newNode->prob = prob;
    newNode->next = NULL;
    return newNode;
}

// A utility function that creates a graph of V vertices
Graph* createGraph(int V)
{
    Graph* graph = new Graph();
    graph->V = V;

    // Create an array of adjacency lists.
    // Size of array will be V
    graph->array = new AdjList[V];

    // Initialize each adjacency list as empty by making head as NULL
    for (int i = 0; i < V; ++i)
        graph->array[i].head = NULL;

    return graph;
}


void printGraph(Graph* graph)
{
    cout << endl << "Printing graph : " << endl;
    int V = graph->V;

    for (int i = 0; i < V; ++i)
    {
        AdjListNode* pointer = graph->array[i].head;
        cout << i << ":  ";
        while (pointer != NULL)
        {
            cout << "(" << pointer->dest << ", " << pointer->weight << ") ";
            pointer = pointer->next;


        }
        cout << endl;
    }
}


int NO_PARENT = -1;


 //dijkstra
// Structure to represent a min heap node
struct MinHeapNode
{
    int v;
    int dist;
};

// Structure to represent a min heap
struct MinHeap
{
    int size;   // Number of heap nodes present currently
    int capacity; // Capacity of min heap
    int *pos; // This is needed for decreaseKey()
    struct MinHeapNode **array;

    ~MinHeap()
    {
        for (int i = 0; i < size; ++i) {
            delete array[i];
        }
        delete[] array;
        delete[] pos;

    }
};

// A utility function to create a new Min Heap Node
struct MinHeapNode* newMinHeapNode(int v,
                                int dist)
{
    MinHeapNode* minHeapNode = new MinHeapNode();
    minHeapNode->v = v;
    minHeapNode->dist = dist;
    return minHeapNode;
}

// A utility function to create a Min Heap
struct MinHeap* createMinHeap(int capacity)
{
    MinHeap* minHeap = new MinHeap();
    minHeap->pos = new int[capacity];
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = new MinHeapNode*[capacity];
    return minHeap;
}

// A utility function to swap two nodes of min heap.
// Needed for min heapify
void swapMinHeapNode( MinHeapNode** a,
                     MinHeapNode** b)
{
    struct MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

// A standard function to heapify at given idx
// This function also updates position of nodes when they are swapped.
// Position is needed for decreaseKey()
void minHeapify( MinHeap* minHeap, int idx)
{
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < minHeap->size &&
        minHeap->array[left]->dist <
        minHeap->array[smallest]->dist )
    smallest = left;

    if (right < minHeap->size &&
        minHeap->array[right]->dist <
        minHeap->array[smallest]->dist )
    smallest = right;

    if (smallest != idx)
    {
        // The nodes to be swapped in min heap
        MinHeapNode *smallestNode =
            minHeap->array[smallest];
        MinHeapNode *idxNode =
                minHeap->array[idx];

        // Swap positions
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;

        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest],
                        &minHeap->array[idx]);

        minHeapify(minHeap, smallest);
    }
}

// A utility function to check if the given minHeap is empty or not
int isEmpty( MinHeap* minHeap)
{
    return minHeap->size == 0;
}

// Standard function to extract minimum node from heap
struct MinHeapNode* extractMin(struct MinHeap*
                                minHeap)
{
    if (isEmpty(minHeap))
        return NULL;

    // Store the root node
    struct MinHeapNode* root =
                minHeap->array[0];

    // Replace root node with last node
    struct MinHeapNode* lastNode =
        minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    // Update position of last node
    minHeap->pos[root->v] = minHeap->size-1;
    minHeap->pos[lastNode->v] = 0;

    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);

    return root;
}

// Function to decreasekey dist value of a given vertex v. This function
// uses pos[] of min heap to get the current index of node in min heap
void decreaseKey(struct MinHeap* minHeap,
                        int v, int dist)
{
    // Get the index of v in heap array
    int i = minHeap->pos[v];

    // Get the node and update its dist value
    minHeap->array[i]->dist = dist;

    // Travel up while the complete
    // tree is not heapified.
    // This is a O(Logn) loop
    while (i && minHeap->array[i]->dist <
        minHeap->array[(i - 1) / 2]->dist)
    {
        // Swap this node with its parent
        minHeap->pos[minHeap->array[i]->v] =
                                    (i-1)/2;
        minHeap->pos[minHeap->array[
                            (i-1)/2]->v] = i;
        swapMinHeapNode(&minHeap->array[i],
                &minHeap->array[(i - 1) / 2]);

        // move to parent index
        i = (i - 1) / 2;
    }
}

// A utility function to check if a given vertex 'v' is in min heap or not
bool isInMinHeap(struct MinHeap *minHeap, int v)
{
if (minHeap->pos[v] < minHeap->size)
    return true;
return false;
}

// A utility function used to print the solution
void printArr(int dist[], int n)
{
    printf("Vertex Distance from Source\n");
    for (int i = 0; i < n; ++i)
        printf("%d \t\t %d\n", i, dist[i]);
}

// The main function that calculates
// distances of shortest paths from src to all
// vertices. It is a O(ELogV) function
vector<int> dijkstra(Graph* graph, int src)
{
    
    // Get the number of vertices in graph
    int V = graph->V;

    // dist values used to pick minimum weight edge in cut
    int dist[V];    

    // minHeap represents set E
    struct MinHeap* minHeap = createMinHeap(V);

    // Initialize min heap with all vertices. dist value of all vertices
    for (int v = 0; v < V; ++v)
    {
        dist[v] = 9999999;
        minHeap->array[v] = newMinHeapNode(v,dist[v]);
        minHeap->pos[v] = v;
    }

    // Make dist value of src vertex as 0 so that it is extracted first
    minHeap->array[src] = newMinHeapNode(src, dist[src]);
    minHeap->pos[src] = src;
    dist[src] = 0;
    decreaseKey(minHeap, src, dist[src]);

    // Initially size of min heap is equal to V
    minHeap->size = V;

    // In the following loop, min heap contains all nodes
    // whose shortest distance is not yet finalized.
    while (!isEmpty(minHeap))
    {
        // Extract the vertex with
        // minimum distance value
        MinHeapNode* minHeapNode =  extractMin(minHeap);
    
        // Store the extracted vertex number
        int u = minHeapNode->v;

        // Traverse through all adjacent
        // vertices of u (the extracted
        // vertex) and update
        // their distance values
        AdjListNode* pCrawl = graph->array[u].head;
        while (pCrawl != NULL)
        {
            int v = pCrawl->dest;

            // If shortest distance to v is
            // not finalized yet, and distance to v
            // through u is less than its
            // previously calculated distance
            if (isInMinHeap(minHeap, v) &&
                    dist[u] != 999999 &&
            pCrawl->weight + dist[u] < dist[v])
            {
                dist[v] = dist[u] + pCrawl->weight;

                // update distance
                // value in min heap also
                decreaseKey(minHeap, v, dist[v]);
            }
            pCrawl = pCrawl->next;
        }

        delete minHeapNode;
        delete pCrawl;
    }



    // print the calculated shortest distances
    vector<int> distVector;
    for (int i = 0; i < V; i++) {
        distVector.push_back(dist[i]);
    }

    delete minHeap;


    return distVector;
}