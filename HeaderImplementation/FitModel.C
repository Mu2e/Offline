#include "TMath.h"
#include <iostream>
#include "config.hh"

using namespace std;

// Switch Double_t to Double_t


// adc fit functions namespace


// THIS CODE WILL NOT WORK SINCE configStruct is not defined.
configStruct initParams;

// Par0 - Q
// The dynamic pedestal is of the form Q e^(-t / tau). 
//Note that normalized it is Q / tau rather than Q 
Float_t dynamicPedestal(Double_t *x, Double_t *par)
{
  return par[0] * exp(-x[0] / initParams._shapingTime);
}

Float_t fixedTruncation(Float_t currentFunctionValue)
{
  return TMath::Min(currentFunctionValue,(Float_t)initParams._truncationLevel);
}


    // Shaping power set to 1
    // MAYBE GET RID OF PAR
Float_t unConvolvedSinglePeak(Double_t *x, Double_t *par)
{
		// Initial return value
		Float_t returnValue = 0.0;
		
    // Set x value
		Float_t xValue = x[0];
  
		if (xValue > 0.0)
		{
  	 returnValue = xValue*pow(initParams._shapingTime,-2)
                    *exp(-xValue/initParams._shapingTime);
		}
		return returnValue;
}
