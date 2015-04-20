#ifndef PeakFitRootBase_hh
#define PeakFitRootBase_hh

#include "TrkChargeReco/inc/PeakFitBase.hh"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TrkChargeReco/inc/ParamStructs.hh"
#include "TrkChargeReco/inc/FitModelRoot.hh"

namespace mu2e {

  namespace TrkChargeReco {

    // Class providing ROOT specific functions for PeakFit classes
    // Contains purely virtual process method (inherited from PeakFitBase)
    class PeakFitRootBase : public PeakFitBase{
      public:

	// PeakFitBaseRoot normal constructor with ConfigStruct initilization parameters
	PeakFitRootBase(const ConfigStruct &initParams) : PeakFitBase(initParams){};

	// Computes initial "guesses" for the resultant hit data
	// This data is passed to all derived classes of PeakFitBaseRoot to compute initial parameters
	// Calls findPeaks and addEarlyPeak
	void initialPeakGuess(const adcWaveform adcData, resultantHitData &initialGuess);

      protected:

	// Fits a model function to a waveform
	void fitModel2NormalizedWaveform(TF1 &fitModel, TGraphErrors &fitData, const Double_t *initialParameters, Double_t *fitParameters);

	// Converts adcWaveform object to TGraphErrors object for easier manipulation in ROOT
	void adcWaveform2TGraphErrors(const adcWaveform adcData, TGraphErrors &fitData);

	TGraphErrors _fitData;
	TF1 _fitModel;

      private:

	// Performs explicit peak search on adc waveform data
	void findPeaks(const TGraphErrors &gr, const ConfigStruct &initParams, resultantHitData &initialGuess, const double sigma = 3.0);

	// This function searches for another peak in the waveform data by subtracting out a dynamic pedestal 
	// from the adc waveform and finding the maximum adc value in the "subtracted data".
	// This function is applied when no peak is found in the explicit peak search (findPeaks).
	void addEarlyPeak(const TGraphErrors &gr, resultantHitData &initialGuess);


    };
  }
}
#endif
