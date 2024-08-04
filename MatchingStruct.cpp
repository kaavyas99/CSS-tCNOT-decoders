#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <set>
#include "BasicFunctions.h"
#include "DistanceFunctions.hpp"
#include "MatchingStruct.h"

/***
 * The idea here is that we create an adjacency list for each surface code (incl time dimension).
 * Note that this graph is accessible at the TwoSCs level not at the SC level (because our SC is a single-time snapshot).
 * At each gate we update the adjacency list using p = p(1-p_new) + (1-p)*pnew. 
 * Finally after all gates are done we convert the adjacency list from probs to weights using w = -ln(p/1-p)
 * 
 * ***/


MatchingStruct createMatchingGraph(int d, double p, double pErase)
{
    MatchingStruct MS;

    MS.d = d;
    MS.p = p;
    MS.pErase = pErase;

    MS.NoOfGraphAncPerLayer = (d+1)*(d+1)/2;
    MS.NoOfAnc = (d+1)*(d+1);
    MS.NoOfLayers = 2*MS.d+1;
    MS.NoOfGraphNodes = MS.NoOfLayers * MS.NoOfGraphAncPerLayer + 2; //add two boundary nodes

    MS.XMatchingGraph = createGraph(MS.NoOfGraphNodes); //two boundaries on this subgraph, front and back 
    MS.ZMatchingGraph = createGraph(MS.NoOfGraphNodes); //left and right 

    return MS;

}

int probsToWeights(double p)
{
    if (p!=0) return int(-1000*log(p/(1-p)));
    else return 99999991;
}

void MatchingStruct::convertProbstoWeights()
{
    //iterate through all elements of X and Z adjacency lists and convert probs to weights

    int V  = XMatchingGraph->V;

    for (int i = 0; i < V; ++i)
    {
        AdjListNode* pCrawl = XMatchingGraph->array[i].head;

        while (pCrawl != NULL)
        {

            double p = pCrawl->prob;
            pCrawl->weight = probsToWeights(p);

            pCrawl = pCrawl->next;
        }
    }


    int V2  = ZMatchingGraph->V;

    for (int i = 0; i < V2; ++i)
    {
        AdjListNode* pCrawl2 = ZMatchingGraph->array[i].head;

        while (pCrawl2 != NULL)
        {

            double p = pCrawl2->prob;
            pCrawl2->weight = probsToWeights(p);

            pCrawl2 = pCrawl2->next;
        }
    } 

}



void baseUpdateEdgeProbs(Graph* graph, int srcStab, int destStab, double p)
{
    if (p > 0) //only worth updating if probability is added to
    {

        int src = srcStab/2;
        int dest = destStab/2;

        //do src - > dest
        AdjListNode* pCrawl = graph->array[src].head;
        bool existsSD = false;

        while (pCrawl != NULL and !existsSD)
        {
            int v = pCrawl->dest;

            if (v == dest)
            {
                existsSD = true;

               double prevp = pCrawl->prob;
               pCrawl->prob = prevp * (1-p) + (1-prevp) * p;

            }
            pCrawl = pCrawl->next;
        }


        //now do dest -> src

        AdjListNode* pCrawl2 = graph->array[dest].head;
        bool existsDS = false;

        while (pCrawl2 != NULL and !existsDS)
        {
            int v = pCrawl2->dest;

            if (v == src)
            {
                existsDS = true;

               double prevp = pCrawl2->prob;
               pCrawl2->prob = prevp * (1-p) + (1-prevp) * p;

            }
            pCrawl2 = pCrawl2->next;
        }


        if (!existsSD)
        {
            baseAddEdge(graph, src, dest, p);
        }

        if (!existsDS and existsSD)
        {
            cout << "Something wrong " << endl;
        }

    }

}

void InsertZerosOfGraph1IntoGraph2(Graph* graph1, Graph* graph2, int startingNode)
{
    for (int i = startingNode; i < graph1->V; ++i) 
    {
        AdjListNode* pCrawl1 = graph1->array[i].head;

        while (pCrawl1 != NULL) 
        {
            // Check if the edge probability in graph1 is zero
            if (pCrawl1->prob == 0.5 and pCrawl1->dest >=startingNode) 
            {
                AdjListNode* pCrawl2 = graph2->array[i].head; // Find the corresponding edge in graph2
                bool foundEdge = false;

                while (pCrawl2 != NULL) 
                {
                    if (pCrawl2->dest == pCrawl1->dest) { foundEdge = true; break; }
                    pCrawl2 = pCrawl2->next;
                }

                // If the edge doesn't exist in graph2, initialize and set it to zero
                if (!foundEdge) 
                {
                    AdjListNode* pCrawl2 = newAdjListNode(pCrawl1->dest, 0.5);
                    pCrawl2->next = graph2->array[i].head;
                    graph2->array[i].head = pCrawl2;

                }
                else { pCrawl2->prob = 0.5; }
            }

            pCrawl1 = pCrawl1->next;
        }
    }
}

void baseAddEdge(Graph* graph, int src, int dest, double prob)
{
	    // Add an edge from src to dest.
        // A new node is added to the adjacency list of src at the beginning
        struct AdjListNode* newNode = newAdjListNode(dest, prob);
        newNode->next = graph->array[src].head;
        graph->array[src].head = newNode;

        // Since graph is undirected, add an edge from dest to src also
        newNode = newAdjListNode(src, prob);
        newNode->next = graph->array[dest].head;
        graph->array[dest].head = newNode;
}



