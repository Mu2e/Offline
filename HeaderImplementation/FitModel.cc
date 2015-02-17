#include "TMath.h"
#include "ParamStructs.hh"
#include "FitModel.hh"
#include <iostream>


//TODO : FINISH TRUNCATION FIT FUNCTIONS

namespace FitModel
{

  Float_t dynamicPedestal(const Double_t t, const DynamicPedestalParamStruct &fitParams, const ConfigStruct &initParams)
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

  Float_t singlePeakWithConstantPedestal(const Double_t t, const SinglePeakWithConstantPedestalParamStruct &fitParams, const ConfigStruct &initParams)
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
  Float_t singlePeakWithDynamicPedestal(const Double_t t, const SinglePeakWithDynamicPedestalParamStruct &fitParams, const ConfigStruct &initParams)
  {
    SinglePeakParamStruct singlePeakFitParams(fitParams._shiftedTime, fitParams._scalingFactor, fitParams._sigma);

    DynamicPedestalParamStruct dynamicPedestalFitParams(fitParams._Q);

    return (Float_t) singlePeak(t, singlePeakFitParams, initParams) 
                + dynamicPedestal(t, dynamicPedestalFitParams, initParams);
  }

  // PUT THE DOUBLE_t PEAK FUNCTIONS HERE

      // Par0 - shift in X 1st peak
      // Par1 - scalingFactor 1st peak
      // Par2 - shift in 2nd peak minus shift in 1st peak
      // Par3 - scaling factor 2nd peak
  Float_t doublePeak(const Double_t t, const DoublePeakParamStruct &fitParams, const ConfigStruct &initParams)
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
  Float_t doublePeakWithConstantPedestal(const Double_t t, const DoublePeakWithConstantPedestalParamStruct &fitParams, const ConfigStruct &initParams)
  {
      DoublePeakParamStruct doublePeakFitPar(fitParams._shiftedTimeFirstPeak, fitParams._scalingFactorFirstPeak,
                                      fitParams._shiftedTimeSecondPeak, fitParams._scalingFactorSecondPeak);

      return doublePeak(t, doublePeakFitPar, initParams) + fitParams._verticalShift;
  }

  // Par0 - shift in X 1st peak
  // Par1 - scalingFactor 1st peak
  // Par2 - Q
  // Par3 - shift in 2nd peak minus shift in 1st peak
  // Par4 - scaling factor 2nd peak
  Float_t doublePeakWithDynamicPedestal(const Double_t t, const DoublePeakWithDynamicPedestalParamStruct &fitParams, const ConfigStruct &initParams)
  {
    DoublePeakParamStruct doublePeakFitPar(fitParams._shiftedTimeFirstPeak, fitParams._scalingFactorFirstPeak, 
                                      fitParams._shiftedTimeSecondPeak, fitParams._scalingFactorSecondPeak);

    DynamicPedestalParamStruct dynamicPedestalFitPar(fitParams._Q);

    return doublePeak(t, doublePeakFitPar, initParams) 
        + dynamicPedestal(t, dynamicPedestalFitPar, initParams);
  }

  /**Float_t fitModel2ADC(const Double_t voltageValue, const ConfigStruct &initParams)
  {
    return  fixedTruncation(initParams._bits2scalingFactor * voltageValue);
  }**/



  // Apply truncation to necessary fit models

  Float_t dynamicPedestalTrunc(const Double_t t, const DynamicPedestalParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(dynamicPedestal(t,fitParams,initParams), initParams);
  }

  Float_t singlePeakTrunc(const Double_t t, const SinglePeakParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(singlePeak(t,fitParams,initParams), initParams);
  }

  Float_t singlePeakWithConstantPedestalTrunc(const Double_t t, const SinglePeakWithConstantPedestalParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(singlePeakWithConstantPedestal(t,fitParams,initParams), initParams);
  }

  Float_t singlePeakWithDynamicPedestalTrunc(const Double_t t, const SinglePeakWithDynamicPedestalParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(singlePeakWithDynamicPedestal(t,fitParams,initParams), initParams);
  }

  Float_t doublePeakTrunc(const Double_t t, const DoublePeakParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(doublePeak(t,fitParams,initParams), initParams);
  }

  Float_t doublePeakWithConstantPedestalTrunc(const Double_t t, const DoublePeakWithConstantPedestalParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(doublePeakWithConstantPedestal(t,fitParams,initParams), initParams);
  }

  Float_t doublePeakWithDynamicPedestalTrunc(const Double_t t, const DoublePeakWithDynamicPedestalParamStruct &fitParams, const ConfigStruct &initParams)
  {
    return fixedTruncation(doublePeakWithDynamicPedestal(t,fitParams,initParams), initParams);
  }



}