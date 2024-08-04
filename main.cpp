#include <iostream>
#include <iomanip>      
#include <ctime>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <sstream>
#include "BasicFunctions.h"
#include "DistanceFunctions.hpp"
#include "MatchingStruct.h"
#include "SurfaceCodeStruct.h"
#include "ExtractSyndromeData.h"
#include "MatchingDecoder.h"
using std::cout;

int main(int argc, char** argv)
{   

    // main 0:operation 1:d_min 2:d_max 3:d_step 4:p_min 5:p_max 6:p_step 7:R_e 8:erasure_form 9:NoSamples 10:rand_seed

    // * Set seed *
    srand(atoi(argv[11]));

    char decoderVariant = argv[1][0]; //d for decoupled decoder, s for single update decoder, o for ordered decoding of ctrl first then target
    bool doCNOT = false; if (decoderVariant != 'd') doCNOT = true; //the other decoders DO NOT work for cNOT performed

    cout << "simulation paramters: (transversal_CNOT=" << doCNOT;
    cout << "), (decoder_variant=" << decoderVariant<< ")" << endl;
    
    for (int sysSize = atoi(argv[2]); sysSize < atoi(argv[3]); sysSize+=atoi(argv[4]))
    {

    // * Initialise simulation parameters *        
    int d = sysSize; // System size, odd integer

    // **** Initialise error parameters ***
    for (double pTotal = atof(argv[5]); pTotal < atof(argv[6]); pTotal+=atof(argv[7]))
    {

    double p = pTotal * (1-atof(argv[8])) ; // physical Pauli error rate
    double pErase = pTotal * atof(argv[8]); // erasure rate
    char erasureModel = argv[9][0]; //form of erasure - u is unbiased, n is biased native, b is biased with BCX
   
   
    int FailCount = 0; // count failures 
    int ctrlXFailCount = 0;
    int ctrlZFailCount = 0;
    int trgtXFailCount = 0;
    int trgtZFailCount = 0;

    int NoOfSamples = atoi(argv[10]); // Vary this parameter accordingly


    for (int samples = 0; samples < NoOfSamples; samples++)
    {
        // *** Initialise system ***
        TwoSCsClass bothSurfaceCodes = initializeRecordArrays(d, p, pErase, erasureModel);

        // *** Perfom simulation, simultaneously building matching graph ***
        bothSurfaceCodes.ftTransversalCNOT(doCNOT);

        //*** get stabilizers from measurement outcomes and weights from probs *** 
        int DecoderEstimate = 0;
        bothSurfaceCodes.PostComputationProcessing(decoderVariant);
        MatchingDecoder SampleDecoder = setupDecoder(bothSurfaceCodes, decoderVariant);

        // ** Decode errors and get parity estimate for logical
        SampleDecoder.getMatchingCorrections();
        DecoderEstimate = SampleDecoder.getDecoderEstimate();


        // ** Check if decoding succeeeded
        int Measurement = bothSurfaceCodes.LogicalMeasurementOutcome();
        if (DecoderEstimate != Measurement) { FailCount++;}
        if (DecoderEstimate%2 != Measurement%2) { ctrlXFailCount++; }
        if ((DecoderEstimate%4)/2 != (Measurement%4)/2) { ctrlZFailCount++; }
        if ((DecoderEstimate%8)/4 != (Measurement%8)/4) { trgtXFailCount++; }
        if (DecoderEstimate/8 != Measurement/8) { trgtZFailCount++; }

    
    }

    cout << d << " " << pTotal << " " << argv[8] << " " << erasureModel << " ";
    cout  << ctrlXFailCount << " "  << ctrlZFailCount << " "  << trgtXFailCount << " "  << trgtZFailCount;
    cout << " " << FailCount << " " << NoOfSamples << endl;

    }
    }

    return 0;
}