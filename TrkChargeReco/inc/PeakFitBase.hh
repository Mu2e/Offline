#ifndef TrkChargeReco_PeakFitBase_hh
#define TrkChargeReco_PeakFitBase_hh

#include "TMath.h"
#include "TrkChargeReco/inc/ConfigStruct.hh"

namespace mu2e {

  namespace TrkChargeReco {

    typedef unsigned int * adcWaveform;

    struct resultantPeakData
    {

      Float_t _peakTime; // time of peak relative to 140.0 ns interval
      Float_t _peakHeight; // in units of bits

      // Default constructor 
      resultantPeakData() : _peakTime(0.0), _peakHeight(0.0){};

      // True constructor
      resultantPeakData(Float_t peakTime, Float_t peakHeight) : _peakTime(peakTime), _peakHeight(peakHeight){};
    };

    // This is object top which will be filled by the process method 
    typedef std::vector<resultantPeakData> resultantHitData;

    //  PeakFitBase
    //  Virtual class providing structure for SinglePeakFit, LXPeakFit, ComboPeakFit, etc. 
    class PeakFitBase {
      public:

	// Fills result using adc waveform data
	virtual void process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result) = 0;

	// Destructor
	virtual ~PeakFitBase(){}

	// PeakFitBase normal constructor with ConfigStruct initilization parameters
	PeakFitBase(const ConfigStruct &initParams) : _initParams(initParams){}
      protected:

	ConfigStruct _initParams; 
    };
  }
}
#endif
