#include "TMath.h"
#include <iostream>
#include "TF1.h"

using namespace std;

// Switch Double_t to Double_t

// ADD FIXED TRUNCATION
 

// These need to eventually be DELETED
const Double_t shapingTime = 25.0;
const Double_t truncationLevel = 1023.0 - 64.0;

// Par0 - Q
// The dynamic pedestal is of the form Q e^(-t / tau). 
//Note that normalized it is Q / tau rather than Q 
Float_t dynamicPedestal(Double_t *x, Double_t *par)
{
  return par[0] * exp(-x[0] / shapingTime);
}

Float_t fixedTruncation(Float_t currentFunctionValue)
{
  return TMath::Min(currentFunctionValue,(Float_t)(1023.0-64.0));
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
  	 returnValue = xValue*pow(shapingTime,-2)
                    *exp(-xValue/shapingTime);
		}
		return returnValue;
	}


		// Note that this is a convolution with a uniform distribution
		//2 Parameters (shaping power set to 1.0)
		//par[0] - sigma
Float_t convolvedSinglePeak(Double_t *x, Double_t *par)
{
  	Float_t returnValue = 0.0;
		
		if (par[0] == 0.0)
    	{
      		Double_t parameters[2] = {1.0,shapingTime}; 
      		returnValue = unConvolvedSinglePeak(x,parameters);
    	}
    	else
    	{
      		const Float_t a = TMath::Max((x[0] + par[0]) / shapingTime,0.0);
      		// Assuming that shaping time is positive and thus b is negative (if t - sigma is)
      		const Float_t b = TMath::Max((x[0] - par[0]) / shapingTime,0.0);
      		returnValue =  (-exp(-a)*(1+a) + exp(-b)*(1+b)) / (2.0 * par[0]);
    	}
    	return returnValue;
  }

//Fitting function for current Function2
//par[0] is shifted time 1st peak
//par[1] is scalingfactor 1st peak
//par[2] is sigma 1st peak
Float_t convolutionSinglePeak(Double_t *x, Double_t *par)
{
  Double_t convolvedSinglePeakX[1] = {x[0] - par[0]};
  Double_t convolvedSinglePeakParams[2] = {par[2]};

  return (Float_t) par[1] * convolvedSinglePeak(convolvedSinglePeakX,convolvedSinglePeakParams);

}
			
//Fitting function for current Function2
//par[0] is shifted time 1st peak
//par[1] is scalingfactor 1st peak
//par[2] is vertical shift 1st peak
//par[3] is sigma 1st peak

Float_t convolutionSinglePeakWithConstantPedestal(Double_t *x, Double_t *par)
{
  Double_t convolutionSinglePeakParams[3] = {par[0],par[1],par[3]};

  return (Float_t) convolutionSinglePeak(x,convolutionSinglePeakParams) + par[2];
}

// This is a truncating fitting function with a dynamical pedestal
// par[0] is shifted time
// par[1] is scaling factor
// par[2] is Q
// par[3] is sigma

// This should inherit from convolutionSinglePeak not FitModelBase
Float_t convolutionSinglePeakWithDynamicPedestal(Double_t *x, Double_t *par)
{
  Double_t convolutionSinglePeakParams[4] = {par[0],par[1],0.0,par[3]};

  Double_t dynamicPedestalParam[1] = {par[2]};

  return (Float_t) convolutionSinglePeak(x,convolutionSinglePeakParams) 
              + dynamicPedestal(x,dynamicPedestalParam);
}

// PUT THE DOUBLE_t PEAK FUNCTIONS HERE

    // Par0 - shift in X 1st peak
    // Par1 - scalingFactor 1st peak
    // Par2 - shift in 2nd peak minus shift in 1st peak
    // Par3 - scaling factor 2nd peak
Float_t doublePeak(Double_t *x, Double_t *par)
{
      // Why can't these be constant arrays???
      Double_t firstPeakPar[3] = {par[0], par[1], 0.0};
      Double_t secondPeakPar[3] = {par[2] + par[0], par[3], 0.0};

      return convolutionSinglePeak(x, firstPeakPar) 
           + convolutionSinglePeak(x, secondPeakPar);
}
 
  // Par0 - shift in X 1st peak
  // Par1 - scalingFactor 1st peak
  // Par2 - vertical shift
  // Par3 - shift in 2nd peak minus shift in 1st peak
  // Par4 - scaling factor 2nd peak
Float_t doublePeakWithConstantPedestal(Double_t *x, Double_t *par)
{
    Double_t doublePeakParams[4] = {par[0],par[1],par[3],par[4]};
    return doublePeak(x, doublePeakParams) + par[2];
}

// Par0 - shift in X 1st peak
// Par1 - scalingFactor 1st peak
// Par2 - Q
// Par3 - shift in 2nd peak minus shift in 1st peak
// Par4 - scaling factor 2nd peak
Float_t doublePeakWithDynamicPedestal(Double_t *x, Double_t *par)
{
  Double_t doublePeakParams[4] = {par[0],par[1],par[3],par[4]};
  Double_t dynamicPedestalParam[1] = {par[2]};

  return doublePeak(x,doublePeakParams) 
      + dynamicPedestal(x,dynamicPedestalParam);
}
