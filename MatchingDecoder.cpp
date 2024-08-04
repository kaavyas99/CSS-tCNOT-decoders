
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <algorithm>
#include "BasicFunctions.h"
#include "DistanceFunctions.hpp"
#include "MatchingStruct.h"
#include "SurfaceCodeStruct.h"
#include "ExtractSyndromeData.h"
#include "PerfectMatching.h"
#include "GeomPerfectMatching.h"
#include "MatchingDecoder.h"


MatchingDecoder setupDecoder(TwoSCsClass& SurfaceCodes, char decoderVariant)
{
    MatchingDecoder MD;

    MD.d = SurfaceCodes.d;
    MD.p = SurfaceCodes.p;
    MD.pErase = SurfaceCodes.pErase;

    MD.decoderVariant = decoderVariant; // d,s, or o

    MD.NoOfLayers = SurfaceCodes.NoOfLayers;
    MD.OneSCNoOfAncillas = SurfaceCodes.OneSCNoOfAncillas;

    MD.oneGraphNoOfRegNodes = SurfaceCodes.OneSCNoOfAncillas * (SurfaceCodes.NoOfLayers - 1) /2;
    MD.oneSCNoOfStabilizers = SurfaceCodes.OneSCNoOfAncillas * (SurfaceCodes.NoOfLayers - 1);

    MD.ctrlStabRecord = new bool[SurfaceCodes.OneSCNoOfAncillas  * (SurfaceCodes.NoOfLayers-1)];
    MD.trgtStabRecord = new bool[SurfaceCodes.OneSCNoOfAncillas  * (SurfaceCodes.NoOfLayers-1)];

    for (int l=0; l < SurfaceCodes.NoOfLayers-1;l++)
    {
        for (int a=0; a< SurfaceCodes.OneSCNoOfAncillas; a++)
        {
            int stab = a + l*SurfaceCodes.OneSCNoOfAncillas;
            bool isXStab = SurfaceCodes.controlSC.isXStab(a);

            MD.ctrlStabRecord[stab] = SurfaceCodes.ctrlStabRecord[stab];
            MD.trgtStabRecord[stab] = SurfaceCodes.trgtStabRecord[stab];

            if (decoderVariant == 's' or decoderVariant=='d') //update defect collections for specific decoders
            {
                if (SurfaceCodes.ctrlStabRecord[stab]) {MD.graphAnyons[!isXStab].push_back(stab);}
                if (SurfaceCodes.trgtStabRecord[stab]) {MD.graphAnyons[!isXStab+2].push_back(stab);}                 
            }
            else if (decoderVariant == 'o')
            {
                if (SurfaceCodes.ctrlStabRecord[stab] and !isXStab) {MD.graphAnyons[!isXStab].push_back(stab);}
                if (SurfaceCodes.trgtStabRecord[stab] and isXStab) {MD.graphAnyons[!isXStab+2].push_back(stab);}                 
            }

        }

    }

    //saving the different matching problems in an array
    MD.graphWeights[0] = SurfaceCodes.controlSC.MatchingGraph.XMatchingGraph; 
    MD.graphWeights[1] = SurfaceCodes.controlSC.MatchingGraph.ZMatchingGraph;
    MD.graphWeights[2] = SurfaceCodes.targetSC.MatchingGraph.XMatchingGraph;
    MD.graphWeights[3] = SurfaceCodes.targetSC.MatchingGraph.ZMatchingGraph;

    return MD;

}



vector<int> MatchingDecoder::getStabToNodeWeightArray(int Vertex1Location, int WhichMatching) 
{
    //returns the distances to all other nodes from a given node representing a stabilizer
    vector<int> internalDistances;

    if (Vertex1Location < oneSCNoOfStabilizers) //actual nodes
    {

        int Vertex1Node = (Vertex1Location)/2;
        internalDistances =  dijkstra(graphWeights[WhichMatching], Vertex1Node);  

        vector<int> bdryDistances(internalDistances.end() - 2, internalDistances.end());
        internalDistances.pop_back();  internalDistances.pop_back();  

        for (int i = 0; i < oneGraphNoOfRegNodes; i++) { internalDistances.push_back(9999997); }

        int boundarnode = Vertex1Node + oneGraphNoOfRegNodes;

        internalDistances[boundarnode] = min(bdryDistances[0], bdryDistances[1]);  
        if (bdryDistances[0] < bdryDistances[1])
        {
            listOfNodesConnectingToLowerBdry[WhichMatching].push_back(Vertex1Location);
        }        

    }
    else // boundary nodes
    {
        for (int i = 0; i < oneGraphNoOfRegNodes; i++) { internalDistances.push_back(9999992); }
        //fine to do above because this is never saved for boundary nodes
        for (int i = 0; i < oneGraphNoOfRegNodes; i++) { internalDistances.push_back(0); }
     }  

    return internalDistances;

}

int MatchingDecoder::weightArrayToGraphIndex(int stab)
{
    if (stab < oneSCNoOfStabilizers) return stab/2;
    else
    {
        return oneGraphNoOfRegNodes + (stab - oneSCNoOfStabilizers)/2; //bdry
    }

}

void MatchingDecoder::getMatchingCorrections()
{
    if (decoderVariant=='s' or decoderVariant=='d')
    {
        for (int MatchingNo = 0; MatchingNo < 4; MatchingNo++) { doMatching(MatchingNo);}
    }
    else if (decoderVariant == 'o')
    {
        doMatching(1); doMatching(2);
    
        //use decoder pairs to update decoder pairs in  other two matching graphs
        FilterandProcessSubgraphs();

        //do final matchings
        doMatching(0); doMatching(3);

    }
}


// Decoder function
void MatchingDecoder::doMatching(int MatchingNo)
{

    struct PerfectMatching::Options options;

    int OrigAnyons = graphAnyons[MatchingNo].size();
    int TotalAnyons = 2 * OrigAnyons; // Each anyon can map to a unique boundary node
    int Locations[TotalAnyons];

    int v = 0;   
    for(int j = 0; j < OrigAnyons; j++)
    {
        Locations[v] = graphAnyons[MatchingNo][j];
        Locations[v + OrigAnyons] = graphAnyons[MatchingNo][j] + oneSCNoOfStabilizers;
        v++;   
    }

    int PairElement1[OrigAnyons];
    int PairElement2[OrigAnyons];


    // What follows can be regarded as a black box where minimum weight perfect matching is conducted
    // ====================================================

    int node_num = TotalAnyons;
    int edge_num = TotalAnyons*(TotalAnyons - 1) / 2;

    int edges[2*edge_num];
    int weights[edge_num];

    for(int j = 0; j < edge_num; j++)
    {
        edges[2*j] = 0;
        edges[2*j+1] = 0;
        weights[j] = 0;
    }

    int CurrentEdge = 0;


    for(int j = 0; j < TotalAnyons; j++)
    {

        int Vertex1 = j;
        int Vertex1Location = Locations[Vertex1];
    
        int Connections = TotalAnyons - 1 - j;

        vector<int> Vertex1Distances = getStabToNodeWeightArray(Vertex1Location, MatchingNo);
        
    
        for(int k = 0; k < Connections; k++)
        {
            
            int Vertex2 = j + k + 1;
            int Vertex2Location = Locations[Vertex2];
        
            edges[2*CurrentEdge] = Vertex1;
            edges[2*CurrentEdge+1] = Vertex2;

        
            weights[CurrentEdge] = Vertex1Distances[weightArrayToGraphIndex(Vertex2Location)];
           CurrentEdge++;
        }
    }


    PerfectMatching *pm = new PerfectMatching(node_num, edge_num);
    for (int e=0; e<edge_num; e++)
    {
        pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]);
    }

    pm->options = options;
    pm->Solve();

    int SolutionNumber = 0;

    for (int i=0; i<node_num; i++)
    {
        int j = pm->GetMatch(i);
        if (i < j)
        {
            PairElement1[SolutionNumber] = Locations[i];
            PairElement2[SolutionNumber] = Locations[j];

            DecoderPairs[MatchingNo].push_back(Locations[i]);
            DecoderPairs[MatchingNo].push_back(Locations[j]);
        }
    }

    delete pm;
    
}

void MatchingDecoder::FilterandProcessSubgraphs()
{
    //only called in case of ordered decoding: 
    //Firstly find anyon pairs that are strictly below tCNOT level and flip boundary
    //these definitely flip logical and need to be used to flip decoder estimate
    for (int MatchingNo = 1; MatchingNo <3; MatchingNo++)
    {
        for (int j =0;  j < DecoderPairs[MatchingNo].size(); j+=2)
        {
            int Loc1 = DecoderPairs[MatchingNo][j];
            int Loc2 = DecoderPairs[MatchingNo][j+1]; 

            if (Loc1 < d*OneSCNoOfAncillas) //Loc1 is usually the lower numbered stab; we take this pair regardless and assume it is passed to the dependent subgraph
            {
                int dplus1Loc1  = d*OneSCNoOfAncillas + (Loc1 % OneSCNoOfAncillas);
                int dplus1Loc2 = -1;
                if (Loc2 < oneSCNoOfStabilizers) { dplus1Loc2= d*OneSCNoOfAncillas + (Loc2 % OneSCNoOfAncillas);} //if does not match to boundary
                
                if (MatchingNo == 1 ) {
                    trgtStabRecord[dplus1Loc1] ^= 1; 
                    if (dplus1Loc2 != -1) {trgtStabRecord[dplus1Loc2] ^= 1;}
                    } 
                if (MatchingNo == 2 ) {
                    ctrlStabRecord[dplus1Loc1] ^= 1; 
                    if (dplus1Loc2 != -1)  {ctrlStabRecord[dplus1Loc2] ^= 1;}
                }

            }
        }
    }


    //collect final anyon pairs for these two affected matchings
    for (int l=0; l < NoOfLayers-1;l++)
    {
        for (int a=0; a< OneSCNoOfAncillas; a++)
        {
            int stab = l*OneSCNoOfAncillas+a;
            bool isXStab = ((a / (d+1) + a % (d+1)) % 2 == 0);

            if (ctrlStabRecord[stab] and isXStab) {graphAnyons[!isXStab].push_back(stab);}
            if (trgtStabRecord[stab] and !isXStab) {graphAnyons[!isXStab+2].push_back(stab);}                 
        }
    }
}


int MatchingDecoder::getDecoderEstimate()
{
    int DecoderEstimate = 0;
    bool saveOrderedM1Updates = 0;
    bool saveOrderedM2Updates = 0;

    for (int MatchingNo = 0; MatchingNo < 4; MatchingNo++)
    {
        bool lowerBdryFlipped = 0;

        for (int j =0;  j < DecoderPairs[MatchingNo].size(); j+=2)
        {
            int Loc1 = DecoderPairs[MatchingNo][j];
            int Loc2 = DecoderPairs[MatchingNo][j+1];

            if (min(Loc1, Loc2) < oneSCNoOfStabilizers and max(Loc1,Loc2) >= oneSCNoOfStabilizers) //pair connects to bdry
            {
                if(std::find(listOfNodesConnectingToLowerBdry[MatchingNo].begin(), listOfNodesConnectingToLowerBdry[MatchingNo].end(), min(Loc1,Loc2)) != listOfNodesConnectingToLowerBdry[MatchingNo].end()) 
                {
                    lowerBdryFlipped = !lowerBdryFlipped;

                    if (decoderVariant == 'o')
                    {
                        if (min(Loc1, Loc2) < d*OneSCNoOfAncillas and MatchingNo == 1) {saveOrderedM1Updates ^= 1;}
                        if (min(Loc1, Loc2) < d*OneSCNoOfAncillas and MatchingNo == 2) {saveOrderedM2Updates ^= 1;}
                    }
                }
            }
        }

        DecoderEstimate = DecoderEstimate + lowerBdryFlipped* (1 << MatchingNo);

    }

    if (decoderVariant == 's' or decoderVariant == 'o')
    {
        int ctrlXParity = DecoderEstimate%2;
        int ctrlZParity = (DecoderEstimate%4)/2;
        int trgtXParity = (DecoderEstimate%8)/4;
        int trgtZParity = DecoderEstimate/8;

        if (decoderVariant == 's')
        {
            if (ctrlZParity) trgtZParity ^= 1;
            if (trgtXParity) ctrlXParity ^= 1;
        }
        if (decoderVariant == 'o')
        {
            ctrlXParity = ctrlXParity^saveOrderedM2Updates;
            trgtZParity = trgtZParity^saveOrderedM1Updates;
        }

        return 8*trgtZParity + 4*trgtXParity + 2*ctrlZParity + ctrlXParity;         
    }


    return DecoderEstimate;

}