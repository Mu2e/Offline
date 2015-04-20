#ifndef FitModel_hh
#define FitModel_hh

#include "TMath.h"
#include "TrkChargeReco/inc/ParamStructs.hh"
#include "TrkChargeReco/inc/ConfigStruct.hh"

// Contains functions to be used for fitting by the PeakFit classes

namespace mu2e {

  namespace TrkChargeReco {

    namespace FitModel {

      // Model for early peak approximated as exponential decay
      Float_t earlyPeak(const Double_t t, const EarlyPeakParamStruct &fitParams, const ConfigStruct &initParams);

      // Returns the minimum of currentFunction value and truncationLevel 
      Float_t fixedTruncation(const Float_t currentFunctionValue, const ConfigStruct &initParams);

      // Model for a unconvolved single peak
      // Shaping power set to 1
      Float_t unConvolvedSinglePeak(const Double_t t, const ConfigStruct &initParams);

      // Convolution of single peak with uniform distribution
      //par[0] - sigma
      Float_t convolvedSinglePeak(const Double_t t, const Double_t sigma, const ConfigStruct &initParams);

      //Fitting function for single peak 
      //par[0] - shifted time 1st peak
      //par[1] - scalingfactor 1st peak
      //par[2] - sigma 1st peak
      Float_t singlePeak(const Double_t t, const SinglePeakParamStruct &fitParams, const ConfigStruct &initParams);

      //Fitting function for single peak with a floating pedestal 
      //par[0] - shifted time 1st peak
      //par[1] - scalingfactor 1st peak
      //par[2] - floating pedetsal 1st peak
      //par[3] - sigma 1st peak
      Float_t singlePeakFloatingPedestal(const Double_t t, const SinglePeakFloatingPedestalParamStruct &fitParams, const ConfigStruct &initParams);

      // Fitting function for single peak with early extra peak
      // par[0] - shifted time
      // par[1] - scaling factor
      // par[2] - Q
      // par[3] - sigma
      // This should inherit from convolutionSinglePeak not FitModelBase
      Float_t EXPeak(const Double_t t, const EXPeakParamStruct &fitParams, const ConfigStruct &initParams);

      // Fitting function for single peak with late extra peak`
      // par0 - shift in X 1st peak
      // par1 - scalingFactor 1st peak
      // par2 - timeshift in 2nd peak minus timeshift in 1st peak
      // par3 - scaling factor 2nd peak
      Float_t LXPeak(const Double_t t, const LXPeakParamStruct &fitParams, const ConfigStruct &initParams);

      // Fitting function for single peak with late extra peak and floating pedestal 
      // par0 - shift in X 1st peak
      // par1 - scalingFactor 1st peak
      // par2 - vertical shift
      // par3 - shift in 2nd peak minus shift in 1st peak
      // par4 - scaling factor 2nd peak
      Float_t LXPeakFloatingPedestal(const Double_t t, const LXPeakFloatingPedestalParamStruct &fitParams, const ConfigStruct &initParams);

      // Fitting function for early and late extra peaks 
      // par0 - shift in X 1st peak
      // par1 - scalingFactor 1st peak
      // par2 - Q
      // par3 - shift in 2nd peak minus shift in 1st peak
      // par4 - scaling factor 2nd peak
      Float_t ELXPeak(const Double_t t, const DoublePeakWithEarlyPeakParamStruct &fitParams, const ConfigStruct &initParams);

      // Apply truncation to necessary fit models

      Float_t earlyPeakTrunc(const Double_t t, const EarlyPeakParamStruct &fitParams, const ConfigStruct &initParams);

      Float_t singlePeakTrunc(const Double_t t, const SinglePeakParamStruct &fitParams, const ConfigStruct &initParams);

      Float_t singlePeakFloatingPedestalTrunc(const Double_t t, const SinglePeakFloatingPedestalParamStruct &fitParams, const ConfigStruct &initParams);

      Float_t EXPeakTrunc(const Double_t t, const EXPeakParamStruct &fitParams, const ConfigStruct &initParams);

      Float_t LXPeakTrunc(const Double_t t, const LXPeakParamStruct &fitParams, const ConfigStruct &initParams);

      Float_t LXPeakFloatingPedestalTrunc(const Double_t t, const LXPeakFloatingPedestalParamStruct &fitParams, const ConfigStruct &initParams);

      Float_t ELXPeakTrunc(const Double_t t, const DoublePeakWithEarlyPeakParamStruct &fitParams, const ConfigStruct &initParams);

    }
  }
}
#endif
