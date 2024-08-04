#ifndef MS_FNC_H
#define MS_FNC_H

using  std::cout;

struct MatchingStruct
{
	int d;
	double p;
	double pErase;

	int NoOfGraphAncPerLayer;
	int NoOfAnc;
	int NoOfLayers;
	int NoOfGraphNodes;

	Graph* XMatchingGraph;
	Graph* ZMatchingGraph;

	//Function declarations
	void convertProbstoWeights();
	
};

	void baseUpdateEdgeProbs(Graph* graph, int src, int dest, double p);
	void baseAddEdge(Graph* graph, int src, int dest, double p);
	void InsertZerosOfGraph1IntoGraph2(Graph* graph1, Graph* graph2, int startingNode);


MatchingStruct createMatchingGraph(int d, double p, double pErase);

#endif