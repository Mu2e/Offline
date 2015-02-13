#include "FitModelRoot.hh"
#include "TMath.h"
#include "FitModel.hh"

namespace FitModelRoot
{
	// Par0 - Q
	// The dynamic pedestal is of the form Q e^(-t / tau). 
	//Note that normalized it is Q / tau rather than Q 
	Float_t dynamicPedestal(Double_t *x, Double_t *par)
	{
	  DynamicPedestalParamStruct parStruct(par[0]);	
	  return FitModel::dynamicPedestal(x[0],parStruct, initParams);
	}	

	//Fitting function for current Function2
	//par[0] is shifted time 1st peak
	//par[1] is scalingfactor 1st peak
	//par[2] is sigma 1st peak
	Float_t convolutionSinglePeak(Double_t *x, Double_t *par)
	{
	  SinglePeakParamStruct parStruct(par[0],par[1],par[2]);

	  return FitModel::singlePeak(x[0],parStruct,initParams);

	}
				
	//Fitting function for current Function2
	//par[0] is shifted time 1st peak
	//par[1] is scalingfactor 1st peak
	//par[2] is vertical shift 1st peak
	//par[3] is sigma 1st peak

	Float_t convolutionSinglePeakWithConstantPedestal(Double_t *x, Double_t *par)
	{
	  SinglePeakWithConstantPedestalParamStruct parStruct(par[0],par[1],par[2],par[3]);

	  return FitModel::singlePeakWithConstantPedestal(x[0],parStruct,initParams);
	}

	// This is a truncating fitting function with a dynamical pedestal
	// par[0] is shifted time
	// par[1] is scaling factor
	// par[2] is Q
	// par[3] is sigma

	// This should inherit from convolutionSinglePeak not FitModelBase
	Float_t convolutionSinglePeakWithDynamicPedestal(Double_t *x, Double_t *par)
	{

		SinglePeakWithDynamicPedestalParamStruct parStruct(par[0],par[1],par[2],par[3]);

  		return FitModel::singlePeakWithDynamicPedestal(x[0],parStruct,initParams);
  	}

	// Par0 - shift in X 1st peak
    // Par1 - scalingFactor 1st peak
    // Par2 - shift in 2nd peak minus shift in 1st peak
    // Par3 - scaling factor 2nd peak
	Float_t doublePeak(Double_t *x, Double_t *par)
	{
	    
		DoublePeakParamStruct parStruct(par[0],par[1],par[2],par[3]);

      	return FitModel::doublePeak(x[0],parStruct,initParams);
	}
	 
    // Par0 - shift in X 1st peak
    // Par1 - scalingFactor 1st peak
    // Par2 - vertical shift
    // Par3 - shift in 2nd peak minus shift in 1st peak
    // Par4 - scaling factor 2nd peak
	Float_t doublePeakWithConstantPedestal(Double_t *x, Double_t *par)
	{
		DoublePeakWithConstantPedestalParamStruct parStruct(par[0],par[1],par[2],par[3],par[4]);

      	return FitModel::doublePeakWithConstantPedestal(x[0],parStruct,initParams);
	}

	// Par0 - shift in X 1st peak
	// Par1 - scalingFactor 1st peak
	// Par2 - Q
	// Par3 - shift in 2nd peak minus shift in 1st peak
	// Par4 - scaling factor 2nd peak
	Float_t doublePeakWithDynamicPedestal(Double_t *x, Double_t *par)
	{

		DoublePeakWithDynamicPedestalParamStruct parStruct(par[0],par[1],par[2],par[3],par[4]);

  		return FitModel::doublePeakWithDynamicPedestal(x[0],parStruct,initParams);
	}

}