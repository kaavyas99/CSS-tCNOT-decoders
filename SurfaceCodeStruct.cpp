#include <iostream>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <cstring>
#include "BasicFunctions.h"
#include "DistanceFunctions.hpp"
#include "MatchingStruct.h"
#include "SurfaceCodeStruct.h"


using namespace std;

SurfCodeQubits createSurfaceCode(int d, double p, double pErase, char erasureModel)
{
	SurfCodeQubits SurfaceCode;

	SurfaceCode.d = d;
	SurfaceCode.NoOfQubits = d*d;
	SurfaceCode.NoOfAncillas = (d+1) * (d+1);

	SurfaceCode.p = p;
	SurfaceCode.pErase = pErase;
	SurfaceCode.erasureModel = erasureModel;

    SurfaceCode.QubitsX = new bool[d*d];
    SurfaceCode.QubitsZ = new bool[d*d];

    fill_n( SurfaceCode.QubitsX, d*d, 0 ); // all qubit bool array elements set to false
    fill_n( SurfaceCode.QubitsZ, d*d, 0 ); // all qubit bool array elements set to false

    SurfaceCode.AncillasXErr = new bool[(d+1) * (d+1)];
    SurfaceCode.AncillasZErr = new bool[(d+1) * (d+1)];

    fill_n( SurfaceCode.AncillasXErr, (d+1) * (d+1), 0 ); // all bool array elements set to false
    fill_n( SurfaceCode.AncillasZErr, (d+1) * (d+1), 0 ); // all bool array elements set to false

    SurfaceCode.MatchingGraph = createMatchingGraph(d, p, pErase);

    return SurfaceCode;
}

//function series to get qubit coordinates

int SurfCodeQubits::qubXCoord(int q) { return 2 * (q % d);}
int SurfCodeQubits::qubYCoord(int q) { return 2 * (q / d);}
int SurfCodeQubits::ancXCoord(int a) { return 2 * (a % (d+1)) -1;}
int SurfCodeQubits::ancYCoord(int a) { return 2 * (a / (d+1)) - 1; }

//function series to find ancillas near a qubit q

int SurfCodeQubits::swAnc(int q) { return (q/d) * (d+1) + (q%d);}
int SurfCodeQubits::nwAnc(int q) { return swAnc(q) + (d+1);}
int SurfCodeQubits::seAnc(int q) { return swAnc(q) + 1;}
int SurfCodeQubits::neAnc(int q) { return nwAnc(q) + 1;}


//function series to find qubits near an ancilla a

int SurfCodeQubits::neQub(int a) 
{ 
	if (ancYCoord(a) ==  2*d-1 or ancXCoord(a) == 2*d-1) return -1;
	return (a/(d+1)) * d + (a % (d+1));
}

int SurfCodeQubits::nwQub(int a) 
{ 
	if (ancYCoord(a) ==  2*d-1 or ancXCoord(a) == -1) return -1;
	return (a/(d+1)) * d + (a % (d+1) -1);
}

int SurfCodeQubits::seQub(int a) 
{ 
	if (ancYCoord(a) ==  -1 or ancXCoord(a) == 2*d-1) return -1;
	return (a/(d+1) -1) * d + (a % (d+1) );
}

int SurfCodeQubits::swQub(int a) 
{ 
	if (ancYCoord(a) ==  -1 or ancXCoord(a) == -1) return -1;
	return (a/(d+1) -1) * d + (a % (d+1) -1);
}


bool SurfCodeQubits::isValidAncilla(int a)
{
	if (ancYCoord(a) == -1)
	{
		if (a != 0 and a%2 == 0) {return true;}
		return false;
	}
	if (ancYCoord(a) == 2*d-1)
	{
		if (a != (d+1)*(d+1)-1 and a%2 == 1) {return true;}
		return false;
	}
	if (ancXCoord(a) == -1)
	{
		if (a / (d+1) != d and (a / (d+1)) % 2 == 1) {return true;}
		return false;
	}
	if (ancXCoord(a) == 2*d-1)
	{
		if (a / (d+1) != 0 and (a / (d+1)) % 2 == 0) {return true;}
		return false;
	}

	return true;

}

bool SurfCodeQubits::isXStab(int a)
{
	if ((a / (d+1) + a % (d+1)) % 2 == 0) {return true;}
	return false;
}

bool SurfCodeQubits::isNOrderqub(int q)
{
	//gates are done in refl(N) order, as opposed to refl(Z) order on data qubit
		if ((q / d + q % d) % 2 == 0) {return false;}
	return true;

}

int SurfCodeQubits::getAncillaForGate(int qi, int t)
{
	/****
	 * numbering for X stabilizers 
	 * nwA       neA
	 *    3     1
	 *       q
	 *    4     2
	 * swA       seA
	 * ***/
	int a = -1; //the relevant ancilla to which the gate is done; 
	bool qubisNOrder = isNOrderqub(qi);
	if (t == 1) { a = neAnc(qi); }
	if (t == 2) { if (qubisNOrder) {a = nwAnc(qi); } else {a = seAnc(qi); }}
	if (t == 3) { if (qubisNOrder) {a = seAnc(qi); } else {a = nwAnc(qi); }}
	if (t == 4) { a = swAnc(qi); }

	return a;
}


void SurfCodeQubits::oneStepGates(int t, int MeasRound)
{
	//perform gates, propagating existing errors and potentially introducing new ones
	//simultaneously reweigh matching graph

	for (int q = 0; q < NoOfQubits; q++)
	{
		int a = getAncillaForGate(q,t); //the relevant ancilla to which the gate is done; 

		if (isValidAncilla(a))
		{
			int postGateStats = 0;
			double randErase =  (double)rand() / (double)RAND_MAX;
			bool erased = (randErase < pErase);


			//propagate errors through gates
			if (!isXStab(a)) { CZErrorProp(AncillasXErr[a], AncillasZErr[a], QubitsX[q], QubitsZ[q]); }
			else { CXErrorProp(AncillasXErr[a], AncillasZErr[a], QubitsX[q], QubitsZ[q]);}

			//introduce errors post gates, may be Paulis or erasures
			if (erased)
			{
				PropagateErasureErrors(q,a,t,MeasRound);
			}
			else
			{
				MatchingGraphUpdateProbs(q, a, t, MeasRound, 4*p/15, 'u');
				TwoQubitGatePauliIntroduce(p, AncillasXErr[a], AncillasZErr[a], QubitsX[q], QubitsZ[q]);
			}

		}
	}
}

void SurfCodeQubits::PropagateErasureErrors(int q, int a, int t, int MeasRound)
{
	//erasure flag has gone off. now create errors on qubits and reweigh edges according to error model
	if (erasureModel == 'u')
	{
		MatchingGraphUpdateProbs(q, a, t, MeasRound, 0.5, 'u');
		TwoQubitGatePauliIntroduce(15.0/16.0, AncillasXErr[a], AncillasZErr[a], QubitsX[q], QubitsZ[q]);					
	}
	else //biased erasure 
	{
		double err = (double)rand() / (double)RAND_MAX;
		if (err < 0.5) { AncillasZErr[a] = AncillasZErr[a] ^ 1; }	//ancilla always has possible Z regardless of cnot,cx, native or bcx

		if (erasureModel == 'n' and isXStab(a)) //native gates x stabs -> put x errs on data qubits
		{
			MatchingGraphUpdateProbs(q, a, t, MeasRound, 0.5, 'n');
			if (err < 0.25 or err >= 0.75) { QubitsX[q] = QubitsX[q] ^ 1;}
		}

		else //gates are bias preserving, put Z errs on qubits
		{
			MatchingGraphUpdateProbs(q, a, t, MeasRound, 0.5, 'b');
			if (err < 0.25 or err >= 0.75) { QubitsZ[q] = QubitsZ[q] ^ 1; }
		}
					
	}

}



vector<vector<int> > SurfCodeQubits::PropagateDummyErrors(int q, int a, int t, int trialErr)
{
	/***
	 * This function creates a dummy mini surface code instance for a particular error
	 * It then propagates the error locally through this instance to see which defects it creates
	 * collects the defect locs, and then passes them onto the next function that builds the graph
	 * ***/

	    bool DummyQubXErrs[d * d]; memset(DummyQubXErrs, 0, sizeof(DummyQubXErrs));
	    bool DummyQubZErrs[d * d]; memset(DummyQubZErrs, 0, sizeof(DummyQubZErrs));
	    bool DummyAncXErrs[(d+1)*(d+1)]; memset(DummyAncXErrs, 0, sizeof(DummyAncXErrs));
	    bool DummyAncZErrs[(d+1)*(d+1)]; memset(DummyAncZErrs, 0, sizeof(DummyAncZErrs));

	    bool OneRoundErrRecord[(d+1)*(d+1)]; memset(OneRoundErrRecord, 0, sizeof(OneRoundErrRecord));

		DummyAncZErrs[a] = (trialErr / 3 == 0 and trialErr % 3 != 1); //set ZI,IZ, or ZZ
		DummyQubZErrs[q] = (trialErr / 3 == 0 and trialErr % 3 != 0);

		DummyAncXErrs[a] = (trialErr / 3 == 1 and trialErr % 3 != 1); //set XI,IX, or XX
		DummyQubXErrs[q] = (trialErr / 3 == 1 and trialErr % 3 != 0);

		vector<int> ZSyndromeLocs; vector<int> XSyndromeLocs;
		int offset = 0;

		while (t < 9) //iterate through two layers of gates
		{
			int nextGate = (t+1)%4; 
			if (t % 4 == 3) nextGate = 4;

			if (nextGate == 1) //record measurements after 4th round of gates
			{
				for (int i = 0; i < NoOfAncillas; i++) 
				{
					if (max(abs(ancXCoord(a)-ancXCoord(i)), abs(ancYCoord(a)-ancYCoord(i))) < 5) //don't bother iterating through ancillas too far away
					{
						  if ((DummyAncZErrs[i] and !OneRoundErrRecord[i]) or (!DummyAncZErrs[i] and OneRoundErrRecord[i]))
						  {
							  	if (isXStab(i)) { XSyndromeLocs.push_back(i+offset); }
							  	else { ZSyndromeLocs.push_back(i+offset); }
							  	OneRoundErrRecord[i] = OneRoundErrRecord[i] ^ 1; //record places where syndrome flips
						  }
					}
					DummyAncZErrs[i] = 0; DummyAncXErrs[i] = 0; //reset ancillas
				}
				offset+=NoOfAncillas;
			}


			// propagate errors through gate round
			for (int qi = 0; qi < NoOfQubits; qi++)
			{
				if (max(abs(qubXCoord(q)-qubXCoord(qi)), abs(qubYCoord(q)-qubYCoord(qi))) < 5) //don't bother iterating through qubits too far away
				{
					int ai = getAncillaForGate(qi, nextGate);
					if (isValidAncilla(ai))
					{
						if (!isXStab(ai)) { CZErrorProp(DummyAncXErrs[ai], DummyAncZErrs[ai], DummyQubXErrs[qi], DummyQubZErrs[qi]);} 
						else  {  CXErrorProp(DummyAncXErrs[ai], DummyAncZErrs[ai], DummyQubXErrs[qi], DummyQubZErrs[qi]); }
					}
				}
			}


			t++;
		}

		if (XSyndromeLocs.size() > 2 or ZSyndromeLocs.size() > 2) cout << "TOO BIG" << endl; //these outputs should never be printed

		if (XSyndromeLocs.size() == 1) XSyndromeLocs.push_back(-1);
		if (ZSyndromeLocs.size() == 1) ZSyndromeLocs.push_back(-1);

		vector<vector<int> > BothSyndromeLocs;
		BothSyndromeLocs.push_back(XSyndromeLocs);
		BothSyndromeLocs.push_back(ZSyndromeLocs);
		return BothSyndromeLocs;		
}


void SurfCodeQubits::MatchingGraphUpdateProbs(int q, int a, int t, int MeasRound, double edgeProb, char errSubset)
{
	// try isolated combinations of ZI,IZ,ZZ and XI, IX,XX errors right after the gate between  a-q
	// collect syndrome data and use that to build matching graph

	//qualify err subset again -> 0 is all possible combinations, 1 is from native, {ZI, IX}, 2 is Zs only {ZI, IZ, ZZ}

	int offset = (MeasRound-1)*NoOfAncillas;

	vector<int> errsToIterate;
	switch (errSubset)
	{
		case 'u': 
			errsToIterate.push_back(0); errsToIterate.push_back(1); errsToIterate.push_back(2);
			errsToIterate.push_back(3); errsToIterate.push_back(4); errsToIterate.push_back(5); break;
		case 's': //single ubiased
			errsToIterate.push_back(1); errsToIterate.push_back(4); break; //z and x on data qubit
		case 'n':
			errsToIterate.push_back(0); errsToIterate.push_back(4); break;
		case 'b':
			errsToIterate.push_back(0); errsToIterate.push_back(1); break;
		case 'z':
			errsToIterate.push_back(1); break; //Z on data qubit
		case 'x':
			errsToIterate.push_back(4); break; //X on data qubit
	}


	for (int trialErr = 0; trialErr < errsToIterate.size(); trialErr++)
	{
		vector<vector<int> > BothSyndromeLocs = PropagateDummyErrors(q, a, t, errsToIterate[trialErr]);

		for (int whichGraph = 0; whichGraph < 2; whichGraph++)
		{
			if (BothSyndromeLocs[whichGraph].size() > 1)
			{
				//break down X syndrome to get pairs and add edge
				int s1 = BothSyndromeLocs[whichGraph][0]; int s2 = BothSyndromeLocs[whichGraph][1];
				if (s2 == -1) {s2 = GetBdryPairOfAnc(s1 % NoOfAncillas) + 2*MatchingGraph.NoOfGraphNodes - offset; }
				
				if (whichGraph == 0) baseUpdateEdgeProbs(MatchingGraph.XMatchingGraph, s1+offset, s2+offset, edgeProb);
				if (whichGraph == 1) baseUpdateEdgeProbs(MatchingGraph.ZMatchingGraph, s1+offset, s2+offset, edgeProb);

			}
		}

	}

}


int SurfCodeQubits::GetBdryPairOfAnc(int a)
{
	if (isXStab(a))
	{
		if (ancXCoord(a) == 2*d - 3) return -2;
		if (ancXCoord(a) == 1) return -4;
		else 
		{
			cout << "A non-edge ancilla is connecting to the boundary!" << endl;
			return 3;
		}

	}
	else
	{
		if (ancYCoord(a) == 2*d - 3) return -2;
		if (ancYCoord(a) == 1) return -4;
		else 
		{
			cout << "A non-edge ancilla is connecting to the boundary!" << endl;
			return 3;
		}
	}
}




// Function to measure plaquette operator 
bool SurfCodeQubits::PerfectStabilizerMeasurement(int a)
{

	if (isValidAncilla(a))
	{
		bool outcome = 0;

		int neighbors[4] = {neQub(a), nwQub(a), seQub(a), swQub(a)};

		for (int i=0; i < 4; i++) 
		{
			if (neighbors[i] != -1) //if neighbour is a valid qubit
			{
				if (isXStab(a)) outcome = outcome ^ QubitsZ[neighbors[i]];
				else outcome = outcome ^ QubitsX[neighbors[i]];

			}
		}
		return outcome;

	}
	else return 0;

}