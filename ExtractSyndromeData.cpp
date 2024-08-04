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
#include "SurfaceCodeStruct.h"
#include "ExtractSyndromeData.h"




TwoSCsClass initializeRecordArrays(int d, double p, double pErase, char erasureModel)
{
	TwoSCsClass MeasRecord; //the full record of the operation on two surface codes

	MeasRecord.controlSC = createSurfaceCode(d, p, pErase, erasureModel);
	MeasRecord.targetSC = createSurfaceCode(d, p, pErase, erasureModel);

	MeasRecord.d = d;
	MeasRecord.p = p;
	MeasRecord.pErase = pErase;
	MeasRecord.erasureModel = erasureModel;

	MeasRecord.OneSCNoOfQubits = d*d;
	MeasRecord.OneSCNoOfAncillas = (d+1) * (d+1);

	MeasRecord.NoOfLayers = 2 + 2 * d; //perfect initialization, 2d rounds of stab msmts, and then perfect stab msmt

	MeasRecord.ctrlMeasurementRecord = new bool[MeasRecord.OneSCNoOfAncillas  * MeasRecord.NoOfLayers];
    fill_n( MeasRecord.ctrlMeasurementRecord, MeasRecord.OneSCNoOfAncillas  * MeasRecord.NoOfLayers, 0 ); // all bool array elements set to false

	MeasRecord.trgtMeasurementRecord = new bool[MeasRecord.OneSCNoOfAncillas  * MeasRecord.NoOfLayers];
    fill_n( MeasRecord.trgtMeasurementRecord, MeasRecord.OneSCNoOfAncillas  * MeasRecord.NoOfLayers, 0 ); // all bool array elements set to false

    MeasRecord.ctrlStabRecord = new bool[MeasRecord.OneSCNoOfAncillas  * (MeasRecord.NoOfLayers-1)];
	MeasRecord.trgtStabRecord = new bool[MeasRecord.OneSCNoOfAncillas  * (MeasRecord.NoOfLayers-1)];

    return MeasRecord;
}

void TwoSCsClass::oneRoundofGates(int MeasRound)
{

   //sw gates from ancilla
    controlSC.oneStepGates(1, MeasRound);
    targetSC.oneStepGates(1, MeasRound);

    //use N,Z orientations  for step 2 gates
    controlSC.oneStepGates(2, MeasRound);
    targetSC.oneStepGates(2, MeasRound);


    //use N,Z orientations for step 3 gates
    controlSC.oneStepGates(3, MeasRound);
    targetSC.oneStepGates(3, MeasRound);


    //ne gates from ancilla
    controlSC.oneStepGates(4, MeasRound);
    targetSC.oneStepGates(4, MeasRound);
}


void TwoSCsClass::postRoundMeasurementCollection(int round)
{
	for (int i =0; i < OneSCNoOfAncillas; i++)
	{
		ctrlMeasurementRecord[round*OneSCNoOfAncillas + i] = controlSC.AncillasZErr[i]; //collect ancilla msmts
		trgtMeasurementRecord[round*OneSCNoOfAncillas + i] = targetSC.AncillasZErr[i];

		controlSC.AncillasZErr[i] = 0; //reset ancillas
		targetSC.AncillasZErr[i] = 0;

		controlSC.AncillasXErr[i] = 0; //reset ancillas
		targetSC.AncillasXErr[i] = 0;
	}

}

void TwoSCsClass::FinalPerfectStabilizerUpdate()
{
	for (int i =0; i < OneSCNoOfAncillas; i++)
	{
		ctrlMeasurementRecord[(NoOfLayers-1)*OneSCNoOfAncillas + i] = controlSC.PerfectStabilizerMeasurement(i);
		trgtMeasurementRecord[(NoOfLayers-1)*OneSCNoOfAncillas + i] = targetSC.PerfectStabilizerMeasurement(i);
	}
}


void TwoSCsClass::ftTransversalCNOT(bool doCNOT)
{

	for (int round = 1; round <d+1; round++)
	{
		oneRoundofGates(round);
		postRoundMeasurementCollection(round);
	}


	if (doCNOT) // do tCNOT
	{
		for (int q = 0; q < OneSCNoOfQubits; q++)
		{
			//do transversal CNOT
			CXErrorProp(controlSC.QubitsX[q], controlSC.QubitsZ[q], targetSC.QubitsX[q], targetSC.QubitsZ[q]);
			ApplyandReweightCNOTErrors(q);

		}
	}


	for (int round = d+1; round <2*d+1; round++)
	{
		oneRoundofGates(round);
		postRoundMeasurementCollection(round);
	}

	FinalPerfectStabilizerUpdate();
}


void TwoSCsClass::ApplyandReweightCNOTErrors(int q) // external function for reweighing SCs for tCNOT gates
{
	if ((double)rand() / (double)RAND_MAX >= pErase) //pauli errors on this gate
	{
		TwoQubitGatePauliIntroduce(p, controlSC.QubitsX[q], controlSC.QubitsZ[q], targetSC.QubitsX[q], targetSC.QubitsZ[q]);
		controlSC.MatchingGraphUpdateProbs(q, controlSC.getAncillaForGate(q, 4), 4, d, 8*p/15, 's');
		targetSC.MatchingGraphUpdateProbs(q, targetSC.getAncillaForGate(q, 4), 4, d, 8*p/15, 's');
	}
	else //erasure on this gate
	{
		if (erasureModel == 'u')  //unbiased erasure
		{ 
			TwoQubitGatePauliIntroduce(15.0/16.0, controlSC.QubitsX[q], controlSC.QubitsZ[q], targetSC.QubitsX[q], targetSC.QubitsZ[q]);	
			controlSC.MatchingGraphUpdateProbs(q, controlSC.getAncillaForGate(q, 4), 4, d, 0.5, 's');
			targetSC.MatchingGraphUpdateProbs(q, targetSC.getAncillaForGate(q, 4), 4, d, 0.5, 's');
		}
		else //biased erasure 
		{
			double err = (double)rand() / (double)RAND_MAX;
			if (err < 0.5) { controlSC.QubitsZ[q] = controlSC.QubitsZ[q] ^ 1; }
			controlSC.MatchingGraphUpdateProbs(q, controlSC.getAncillaForGate(q, 4), 4, d, 0.5, 'z');

			if (erasureModel == 'n')  
			{ 
				if (err < 0.25 or err >= 0.75) { targetSC.QubitsX[q] = targetSC.QubitsX[q] ^ 1;} 
				targetSC.MatchingGraphUpdateProbs(q, targetSC.getAncillaForGate(q, 4), 4, d, 0.5, 'x'); //put X erasure on targets
			} 
			else if (err < 0.25 or err >= 0.75) 
			{ 
				targetSC.QubitsZ[q] = targetSC.QubitsZ[q] ^ 1; 
				targetSC.MatchingGraphUpdateProbs(q, targetSC.getAncillaForGate(q, 4), 4, d, 0.5, 'z'); //put Z erasure on targets
			} 
						
		}

	}
	
}

int TwoSCsClass::LogicalMeasurementOutcome() //read out logicals of two SCS
{
	bool ctrlXoutcome = 0;
	bool ctrlZoutcome = 0;
	bool trgtXoutcome = 0;
	bool trgtZoutcome = 0;

	for (int i = 0; i < d; i ++) //we measure along one logical
	{
		ctrlXoutcome = ctrlXoutcome ^ controlSC.QubitsZ[i*d];
		ctrlZoutcome = ctrlZoutcome ^ controlSC.QubitsX[i];
		trgtXoutcome = trgtXoutcome ^ targetSC.QubitsZ[i*d];
		trgtZoutcome = trgtZoutcome ^ targetSC.QubitsX[i];
	}

	return ctrlXoutcome + 2*ctrlZoutcome + 4*trgtXoutcome + 8*trgtZoutcome;
}

void TwoSCsClass::PostComputationProcessing(char decoderVariant)
{

	if (decoderVariant=='s') { posttCNOTErasureWeightUpdater(); } //set any 0 weighted edge on TX graph post tCNOT to 0 in CX graph

	//convert gate probs to weights
	controlSC.MatchingGraph.convertProbstoWeights();
	targetSC.MatchingGraph.convertProbstoWeights();
	
	//printGraph(targetSC.MatchingGraph.ZMatchingGraph);		

	//we have all the measurement outcomes. now we need to convert these to stabilizer outcomes on the relevant matching graphs
	for (int l=0; l <NoOfLayers-1;l++)
	{
		for (int a=0; a<OneSCNoOfAncillas; a++)
		{
			int stab = a + l*OneSCNoOfAncillas;

			bool regCtrlStabLights = ctrlMeasurementRecord[stab] ^ ctrlMeasurementRecord[stab+OneSCNoOfAncillas];
			bool regTrgtStabLights = trgtMeasurementRecord[stab] ^ trgtMeasurementRecord[stab+OneSCNoOfAncillas];

			bool isXStab = controlSC.isXStab(a);
			ctrlStabRecord[stab] = regCtrlStabLights;
			trgtStabRecord[stab] = regTrgtStabLights;

			if (decoderVariant == 's' and l==d) //update depdendent subgraps for dynamic frame
			{
				ctrlStabRecord[stab] = (regCtrlStabLights ^ ( (trgtMeasurementRecord[stab+OneSCNoOfAncillas]) * isXStab));
				trgtStabRecord[stab] = (regTrgtStabLights ^ ( (ctrlMeasurementRecord[stab+OneSCNoOfAncillas] )* !isXStab)); 
			}
			else if (decoderVariant == 's' and l>d) //update depdendent subgraps for dynamic frame
			{
				ctrlStabRecord[stab] = (regCtrlStabLights ^ ( (regTrgtStabLights) * isXStab));
				trgtStabRecord[stab] = (regTrgtStabLights ^ ( (regCtrlStabLights )* !isXStab));
			}

		}

	}

}

void TwoSCsClass::posttCNOTErasureWeightUpdater() //update dependent graph weights
{
    int stabsAftertCNOT = d*OneSCNoOfAncillas/2;

    for(int ZMatchingGraph = 0; ZMatchingGraph < 2; ZMatchingGraph++)
    {
    	Graph* graph2 = controlSC.MatchingGraph.XMatchingGraph;
    	Graph* graph1 = targetSC.MatchingGraph.XMatchingGraph;

    	if (ZMatchingGraph) 
    	{
		   graph1  = controlSC.MatchingGraph.ZMatchingGraph;
		   graph2  = targetSC.MatchingGraph.ZMatchingGraph;
    	}

    	InsertZerosOfGraph1IntoGraph2(graph1, graph2, stabsAftertCNOT);
	}

}
