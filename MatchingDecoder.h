#ifndef MD_FNC_H
#define MD_FNC_H

using namespace std;


struct MatchingDecoder
{
	int d;
	double p;
	double pErase;

	char decoderVariant;

	int oneGraphNoOfRegNodes;
	int oneSCNoOfStabilizers;

	int OneSCNoOfAncillas;
	int NoOfLayers;

	bool* ctrlStabRecord;
	bool* trgtStabRecord;

	vector<int> graphAnyons[4];
	Graph* graphWeights[4];
	vector<int> DecoderPairs[4];

	vector<int> listOfNodesConnectingToLowerBdry[4];


	//function declarations
	int weightArrayToGraphIndex(int stab);
	vector<int> getStabToNodeWeightArray(int Vertex1Location, int WhichMatching);

	void FilterandProcessSubgraphs();

	void getMatchingCorrections();
	void doMatching(int MatchingNo);

	int getDecoderEstimate();

	~MatchingDecoder()
	{
	    // Destructor definitio
	    delete[] ctrlStabRecord;
	    delete[] trgtStabRecord;
	}

	
};


MatchingDecoder setupDecoder(TwoSCsClass& SurfaceCodes, char decoderVariant);



#endif