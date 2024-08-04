#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include<set>

using namespace std;



void CZErrorProp(bool &ControlXErr, bool &ControlZErr, bool &TargetXErr, bool &TargetZErr)
{
    if (ControlXErr)
    {
        TargetZErr = TargetZErr ^ 1;
    }
    if (TargetXErr)
    {
        ControlZErr = ControlZErr ^ 1;
    }

}

void CXErrorProp(bool &ControlXErr, bool &ControlZErr, bool &TargetXErr, bool &TargetZErr)
{
    if (ControlXErr)
    {
        TargetXErr = TargetXErr ^ 1;
    }
    if (TargetZErr)
    {
        ControlZErr = ControlZErr ^ 1;
    }
}

void TwoQubitGatePauliIntroduce(double p, bool &ControlXErr, bool &ControlZErr, bool &TargetXErr, bool &TargetZErr)
{

    double err =  (double)rand() / (double)RAND_MAX;
     
    if (err < p/15)
    {
        TargetXErr = TargetXErr ^ 1;
    }
    else if (err < 2*p/15)
    {
        TargetXErr = TargetXErr ^ 1;
        TargetZErr = TargetZErr ^ 1;
    }
    else if (err < p/5)
    {
        TargetZErr = TargetZErr ^ 1;
    }

    else if (err < 4*p/15)
    {
        ControlXErr = ControlXErr ^ 1;
    }
    else if (err < 5*p/15)
    {
        ControlXErr = ControlXErr ^ 1;
        TargetXErr = TargetXErr ^ 1;
    }
    else if (err < 6*p/15)
    {
        ControlXErr = ControlXErr ^ 1;
        TargetXErr = TargetXErr ^ 1;        
        TargetZErr = TargetZErr ^ 1;
    }
    else if (err < 7*p/15)
    {
        ControlXErr = ControlXErr ^ 1;       
        TargetZErr = TargetZErr ^ 1;
    }

    else if (err < 8*p/15)
    {
        ControlXErr = ControlXErr ^ 1;
        ControlZErr = ControlZErr ^ 1;
    }
    else if (err < 9*p/15)
    {
        ControlXErr = ControlXErr ^ 1;
        ControlZErr = ControlZErr ^ 1;
        TargetXErr = TargetXErr ^ 1;
    }
    else if (err < 10*p/15)
    {
        ControlXErr = ControlXErr ^ 1;
        ControlZErr = ControlZErr ^ 1;
        TargetXErr = TargetXErr ^ 1;        
        TargetZErr = TargetZErr ^ 1;
    }
    else if (err < 11*p/15)
    {
        ControlXErr = ControlXErr ^ 1;
        ControlZErr = ControlZErr ^ 1;     
        TargetZErr = TargetZErr ^ 1;
    }

    else if (err < 12*p/15)
    {
        ControlZErr = ControlZErr ^ 1;
    }
    else if (err < 13*p/15)
    {
        ControlZErr = ControlZErr ^ 1;
        TargetXErr = TargetXErr ^ 1;
    }
    else if (err < 14*p/15)
    {
        ControlZErr = ControlZErr ^ 1;
        TargetXErr = TargetXErr ^ 1;        
        TargetZErr = TargetZErr ^ 1;
    }
    else if (err < p)
    {
        ControlZErr = ControlZErr ^ 1;     
        TargetZErr = TargetZErr ^ 1;
    }
}



