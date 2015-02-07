#ifndef FitModel_hh
#define FitModel_hh

#include "TMath.h"
#include "ParamStructs.hh"
#include "config.hh"

namespace FitModel
{
  Float_t dynamicPedestal(Double_t t, dynamicPedestalParamStruct &fitParams, configStruct &initParams);
  
  Float_t fixedTruncation(Float_t currentFunctionValue, configStruct &initParams);

  // Shaping power set to 1
  Float_t unConvolvedSinglePeak(Double_t t, configStruct &initParams);

	// Note that this is a convolution with a uniform distribution
	//2 Parameters (shaping power set to 1.0)
	//par[0] - sigma
  Float_t convolvedSinglePeak(Double_t t, Double_t sigma, configStruct &initParams);

  //Fitting function for current Function2
  //par[0] is shifted time 1st peak
  //par[1] is scalingfactor 1st peak
  //par[2] is sigma 1st peak
  Float_t singlePeak(Double_t t, singlePeakParamStruct &fitParams, configStruct &initParams);
  			
  //Fitting function for current Function2
  //par[0] is shifted time 1st peak
  //par[1] is scalingfactor 1st peak
  //par[2] is vertical shift 1st peak
  //par[3] is sigma 1st peak

  Float_t singlePeakWithConstantPedestal(Double_t t, singlePeakWithConstantPedestalParamStruct &fitParams, configStruct &initParams);

  // This is a truncating fitting function with a dynamical pedestal
  // par[0] is shifted time
  // par[1] is scaling factor
  // par[2] is Q
  // par[3] is sigma

  // This should inherit from convolutionSinglePeak not FitModelBase
  Float_t singlePeakWithDynamicPedestal(Double_t t, singlePeakWithDynamicPedestalParamStruct &fitParams, configStruct &initParams);

  // PUT THE DOUBLE_t PEAK FUNCTIONS HERE

      // Par0 - shift in X 1st peak
      // Par1 - scalingFactor 1st peak
      // Par2 - shift in 2nd peak minus shift in 1st peak
      // Par3 - scaling factor 2nd peak
  Float_t doublePeak(Double_t t, doublePeakParamStruct &fitParams, configStruct &initParams);
   
    // Par0 - shift in X 1st peak
    // Par1 - scalingFactor 1st peak
    // Par2 - vertical shift
    // Par3 - shift in 2nd peak minus shift in 1st peak
    // Par4 - scaling factor 2nd peak
  Float_t doublePeakWithConstantPedestal(Double_t t, doublePeakWithConstantPedestalParamStruct &fitParams, configStruct &initParams);

  // Par0 - shift in X 1st peak
  // Par1 - scalingFactor 1st peak
  // Par2 - Q
  // Par3 - shift in 2nd peak minus shift in 1st peak
  // Par4 - scaling factor 2nd peak
  Float_t doublePeakWithDynamicPedestal(Double_t t, doublePeakWithDynamicPedestalParamStruct &fitParams, configStruct &initParams);
}
#endif