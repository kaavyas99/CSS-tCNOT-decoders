#ifndef ESD_FNC_H
#define ESD_FNC_H

using namespace std;

struct TwoSCsClass
{
	SurfCodeQubits controlSC;
	SurfCodeQubits targetSC;

	int d; 
	double p;
	double pErase;
	int erasureModel;

	int NoOfQubits;
	int NoOfAncillas;


	int OneSCNoOfQubits;
	int OneSCNoOfAncillas;
	int NoOfLayers;

	bool* ctrlMeasurementRecord;
	bool* trgtMeasurementRecord;

	bool* ctrlStabRecord;
	bool* trgtStabRecord;



	//Function declarations
	void oneRoundofGates(int round);
	void postRoundMeasurementCollection(int round);
	void ftTransversalCNOT(bool doCNOT);
	void ApplyandReweightCNOTErrors(int q);

	void PostComputationProcessing(char decoderVariant);
	void posttCNOTErasureWeightUpdater();

	void FinalPerfectStabilizerUpdate();
	int LogicalMeasurementOutcome();


	~TwoSCsClass()
	{
	    // Destructor definition
	    //cout << ("|");
	    delete[] ctrlMeasurementRecord;
	    delete[] trgtMeasurementRecord;

	    delete[] ctrlStabRecord;
	    delete[] trgtStabRecord;

	    delete[] controlSC.QubitsX;
	    delete[] controlSC.QubitsZ;
	    delete[] targetSC.QubitsX;
	    delete[] targetSC.QubitsZ; 

	   	delete[] controlSC.AncillasXErr;
	    delete[] targetSC.AncillasXErr;
	   	delete[] controlSC.AncillasZErr;
	    delete[] targetSC.AncillasZErr;

	    delete controlSC.MatchingGraph.XMatchingGraph;
	    delete controlSC.MatchingGraph.ZMatchingGraph;

	    delete targetSC.MatchingGraph.XMatchingGraph;
	    delete targetSC.MatchingGraph.ZMatchingGraph;

	}
};

TwoSCsClass initializeRecordArrays(int d, double p, double pErase, char erasureModel);

#endif
