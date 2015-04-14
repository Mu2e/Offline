#include "TMath.h"
#include "TrkChargeReco/inc/ParamStructs.hh"
#include "TrkChargeReco/inc/FitModel.hh"
#include <iostream>


//TODO : FINISH TRUNCATION FIT FUNCTIONS
namespace mu2e {

namespace TrkChargeReco {

namespace FitModel {

  Float_t earlyPeak(const Double_t t, const EarlyPeakParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fitParams._Q * exp(-t / initParams._shapingTime);
  }

  Float_t fixedTruncation(const Float_t currentFunctionValue, const ConfigStruct &initParams)
  {
    Float_t result = currentFunctionValue;
    if (currentFunctionValue > initParams._truncationLevel) 
      {result = initParams._truncationLevel;}
    return result;
  }

  // Shaping power set to 1
  Float_t unConvolvedSinglePeak(const Double_t t, const ConfigStruct &initParams)
  {
  		// Initial return value
  		Float_t returnValue = 0.0;
    
  		if (t > 0.0)
  		{
    	 returnValue = t*pow(initParams._shapingTime,-2)
                      *exp(-t/initParams._shapingTime);
  		}
  		return returnValue;
  	}


  		// Note that this is a convolution with a uniform distribution
  		//2 Parameters (shaping power set to 1.0)
  		//par[0] - sigma
  Float_t convolvedSinglePeak(const Double_t t, const Double_t sigma, const ConfigStruct &initParams)
  {
    	Float_t returnValue = 0.0;
  		
  	if (sigma == 0.0)
      	{
        		returnValue = unConvolvedSinglePeak(t, initParams);
      	}
      	else
      	{
        		const Float_t a = TMath::Max((t + sigma) / initParams._shapingTime,0.0);
        		// Assuming that shaping time is positive and thus b is negative (if t - sigma is)
        		const Float_t b = TMath::Max((t - sigma) / initParams._shapingTime,0.0);
        		returnValue =  (-exp(-a)*(1+a) + exp(-b)*(1+b)) / (2.0 * sigma);
      	}
      	return returnValue;
    }

  //Fitting function for current Function2
  //par[0] is shifted time 1st peak
  //par[1] is scalingfactor 1st peak
  //par[2] is sigma 1st peak
  Float_t singlePeak(const Double_t t, const SinglePeakParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return (Float_t) fitParams._scalingFactor * convolvedSinglePeak(t - fitParams._shiftedTime, fitParams._sigma, initParams) + initParams._defaultPedestal;

  }
  			
  //Fitting function for current Function2
  //par[0] is shifted time 1st peak
  //par[1] is scalingfactor 1st peak
  //par[2] is vertical shift 1st peak
  //par[3] is sigma 1st peak

  Float_t singlePeakFloatingPedestal(const Double_t t, const SinglePeakFloatingPedestalParamStruct &fitParams, const ConfigStruct &initParams)
  {
    SinglePeakParamStruct singlePeakFitParams(fitParams._shiftedTime, fitParams._scalingFactor, fitParams._sigma);
    return (Float_t) singlePeak(t, singlePeakFitParams, initParams) + fitParams._verticalShift;
  }

  // This is a truncating fitting function with a dynamical pedestal
  // par[0] is shifted time
  // par[1] is scaling factor
  // par[2] is Q
  // par[3] is sigma

  // This should inherit from convolutionSinglePeak not FitModelBase
  Float_t EXPeak(const Double_t t, const EXPeakParamStruct &fitParams, const ConfigStruct &initParams)
  {
    SinglePeakParamStruct singlePeakFitParams(fitParams._shiftedTime, fitParams._scalingFactor, fitParams._sigma);

    EarlyPeakParamStruct earlyPeakFitParams(fitParams._Q);

    return (Float_t) singlePeak(t, singlePeakFitParams, initParams) 
                + earlyPeak(t, earlyPeakFitParams, initParams);
  }

  // PUT THE DOUBLE_t PEAK FUNCTIONS HERE

      // Par0 - shift in X 1st peak
      // Par1 - scalingFactor 1st peak
      // Par2 - shift in 2nd peak minus shift in 1st peak
      // Par3 - scaling factor 2nd peak
  Float_t LXPeak(const Double_t t, const LXPeakParamStruct &fitParams, const ConfigStruct &initParams)
  { 
        // Sigma is set to 0 for each peak
        SinglePeakParamStruct firstPeakFitPar(fitParams._shiftedTimeFirstPeak, fitParams._scalingFactorFirstPeak);
       
        // MAYBE CHANGE NAME OF _SHIFTEDTIMESECONDPEAK
        SinglePeakParamStruct secondPeakFitPar(fitParams._shiftedTimeSecondPeak + fitParams._shiftedTimeFirstPeak, 
                                             fitParams._scalingFactorSecondPeak);

        // The default pedestal is subtracted since it is double counted by adding two single peaks
        // together. 
        return singlePeak(t, firstPeakFitPar, initParams) 
             + singlePeak(t, secondPeakFitPar, initParams)
             - initParams._defaultPedestal;
  }
   
    // Par0 - shift in X 1st peak
    // Par1 - scalingFactor 1st peak
    // Par2 - vertical shift
    // Par3 - shift in 2nd peak minus shift in 1st peak
    // Par4 - scaling factor 2nd peak
  Float_t LXPeakFloatingPedestal(const Double_t t, const LXPeakFloatingPedestalParamStruct &fitParams, const ConfigStruct &initParams)
  {
      LXPeakParamStruct LXPeakFitPar(fitParams._shiftedTimeFirstPeak, fitParams._scalingFactorFirstPeak,
                                      fitParams._shiftedTimeSecondPeak, fitParams._scalingFactorSecondPeak);

      return LXPeak(t, LXPeakFitPar, initParams) + fitParams._verticalShift;
  }

  // Par0 - shift in X 1st peak
  // Par1 - scalingFactor 1st peak
  // Par2 - Q
  // Par3 - shift in 2nd peak minus shift in 1st peak
  // Par4 - scaling factor 2nd peak
  Float_t ELXPeak(const Double_t t, const DoublePeakWithEarlyPeakParamStruct &fitParams, const ConfigStruct &initParams)
  {
    LXPeakParamStruct LXPeakFitPar(fitParams._shiftedTimeFirstPeak, fitParams._scalingFactorFirstPeak, 
                                      fitParams._shiftedTimeSecondPeak, fitParams._scalingFactorSecondPeak);

    EarlyPeakParamStruct earlyPeakFitPar(fitParams._Q);

    return LXPeak(t, LXPeakFitPar, initParams) 
        + earlyPeak(t, earlyPeakFitPar, initParams);
  }

  /**Float_t fitModel2ADC(const Double_t voltageValue, const ConfigStruct &initParams)
  {
    return  fixedTruncation(initParams._bits2scalingFactor * voltageValue);
  }**/



  // Apply truncation to necessary fit models

  Float_t earlyPeakTrunc(const Double_t t, const EarlyPeakParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(earlyPeak(t,fitParams,initParams), initParams);
  }

  Float_t singlePeakTrunc(const Double_t t, const SinglePeakParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(singlePeak(t,fitParams,initParams), initParams);
  }

  Float_t singlePeakFloatingPedestalTrunc(const Double_t t, const SinglePeakFloatingPedestalParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(singlePeakFloatingPedestal(t,fitParams,initParams), initParams);
  }

  Float_t EXPeakTrunc(const Double_t t, const EXPeakParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(EXPeak(t,fitParams,initParams), initParams);
  }

  Float_t LXPeakTrunc(const Double_t t, const LXPeakParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(LXPeak(t,fitParams,initParams), initParams);
  }

  Float_t LXPeakFloatingPedestalTrunc(const Double_t t, const LXPeakFloatingPedestalParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(LXPeakFloatingPedestal(t,fitParams,initParams), initParams);
  }

  Float_t ELXPeakTrunc(const Double_t t, const DoublePeakWithEarlyPeakParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(ELXPeak(t,fitParams,initParams), initParams);
  }


}
}
}
