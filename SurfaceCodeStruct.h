#ifndef SC_FNC_H
#define SC_FNC_H

using namespace std;

struct SurfCodeQubits
{
	int d; 
	int NoOfQubits;
	int NoOfAncillas;
	char erasureModel;


	bool* QubitsX;
	bool* QubitsZ;
	bool* AncillasXErr;
	bool* AncillasZErr;

	MatchingStruct MatchingGraph;

	double p;
	double pErase;

	//Function declarations
	int qubXCoord(int q); //find qubit coordinates given qubit index
	int qubYCoord(int q);

	int ancXCoord(int a); //find ancilla coordinates given ancilla index
	int ancYCoord(int a);

	bool isXStab(int a);
	bool isNOrderqub(int q);
	int getAncillaForGate(int qi, int t);


	int swAnc(int q);
	int nwAnc(int q);
	int seAnc(int q);
	int neAnc(int q);

	int neQub(int a); 
	int nwQub(int a); 
	int seQub(int a);
	int swQub(int a);

	bool isValidAncilla(int a);

	void oneStepGates(int t, int MeasRound); //do relevant gate at t=1-4
	vector<vector<int> > PropagateDummyErrors(int q, int a, int t, int trialErr); //helps find syndromes that an error creates
	void PropagateErasureErrors(int q, int a, int t, int MeasRound);
	void MatchingGraphUpdateProbs(int q, int a, int t, int MeasRound, double edgeProb, char errSubset);
	bool PerfectStabilizerMeasurement(int a);
	int GetBdryPairOfAnc(int a);


};

SurfCodeQubits createSurfaceCode(int d, double p, double pErase, char erasureModel);



#endif
