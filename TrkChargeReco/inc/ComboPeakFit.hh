#ifndef ComboPeakFit_hh
#define ComboPeakFit_hh

#include "TrkChargeReco/inc/PeakFitRootBase.hh"

// Contains PeakFit classes 
// To perform fit call the process method
// Note that EXPeakFit and LXPeakFit denote early and late extra peak fits, respectively.
// All peakFit classes currently inherit from the PeakFitRootBase class and hence are required to  have
// an implementation of the process method

namespace mu2e {

  namespace TrkChargeReco {

    // Computes reconstructed energy by summing adc values and subtracting two presamples
    class SumADC : public PeakFitRootBase{
      public:
	SumADC(const ConfigStruct &initParams) : PeakFitRootBase(initParams){};

	virtual void process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result);
    };

    class SinglePeakFit : public PeakFitRootBase{
      public:
	// SinglePeakFit normal constructor with ConfigStruct initilization parameters
	SinglePeakFit(const ConfigStruct &initParams);

	// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
	virtual void process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result);

      protected:

	// Convert parameters from fit to resultantHitData object
	void fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result);
    };

    class SinglePeakFloatingPedestalFit : public PeakFitRootBase{
      public:
	// SinglePeakFit normal constructor with ConfigStruct initilization parameters
	SinglePeakFloatingPedestalFit(const ConfigStruct &initParams);

	// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
	virtual void process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result);

      protected:

	void fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result);
    };

    class EXPeakFit : public PeakFitRootBase{
      public:

	// EXPeakFit normal constructor with ConfigStruct initilization parameters
	EXPeakFit(const ConfigStruct &initParams);

	// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
	virtual void process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result);

      protected:

	// Convert parameters from fit to resultantHitData object
	void fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result);
    };

    class LXPeakFit : public PeakFitRootBase{
      public:
	LXPeakFit(const ConfigStruct &initParams);

	// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
	virtual void process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result);

      protected:

	// Convert parameters from fit to resultantHitData object
	void fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result);
    };

    class LXPeakFloatingPedestalFit : public PeakFitRootBase{
      public:
	LXPeakFloatingPedestalFit(const ConfigStruct &initParams);

	// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
	virtual void process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result);

      protected:

	// Convert parameters from fit to resultantHitData object
	void fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result);
    };

    class ELXPeakFit : public PeakFitRootBase{
      public:
	// ELXPeakFit normal constructor with ConfigStruct initilization parameters
	ELXPeakFit(const ConfigStruct &initParams);

	// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
	virtual void process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result);

      protected:

	// Convert parameters from fit to resultantHitData object
	void fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result);
    };

    class ComboPeakFit : public PeakFitRootBase{
      public:

	// ComboPeakFit normal constructor with ConfigStruct initilization parameters
	ComboPeakFit(const ConfigStruct &initParams) : PeakFitRootBase(initParams){};

	// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
	virtual void process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result); 

    };
  }
}
#endif
