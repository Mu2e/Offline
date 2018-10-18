
#ifndef TrkHitReco_PeakFitRoot_hh
#define TrkHitReco_PeakFitRoot_hh
// base class for controlling peak fits
    //  Virtual class providing structure for extracting charge from ADC waveforms
//
#include "TrkHitReco/inc/PeakFit.hh"
#include "TrkHitReco/inc/PeakFitFunction.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include <string>

class TGraphErrors;


namespace mu2e {

  namespace TrkHitReco {


    class PeakFitRoot : public PeakFit 
    {
      public:
	
	PeakFitRoot(const StrawResponse& srep, const fhicl::ParameterSet& pset);
 	virtual ~PeakFitRoot(){}


	// extract peak information from adc waveform data.  1 waveform generates 1 peak fit.
	// The default implementation simply sums the ADC data after subtracting pedesdal
	virtual void process(TrkTypes::ADCWaveform const& adcData, PeakFitParams & fit) const;
     	
        // Converts adcWaveform object to TGraphErrors object for easier manipulation in ROOT
	void adcWaveform2TGraphErrors(TrkTypes::ADCWaveform const& adcData, TGraphErrors &fitData) const;
      
      protected:
        bool            _truncateADC;     // model ADC truncation
        bool            _floatPedestal;   // float pedestal in fit
        bool            _floatWidth;      // _float width in fit
        bool            _earlyPeak;       // additional peak finding
        bool            _latePeak;        // additional peak finding
	std::string     _fitoptions;      // ROOT fit options
        unsigned        _maxFitIter;      // maximumnumber of fit options
        int             _debug; 
	
	PeakFitFunction _peakfit;
        FitConfig       _config;
    };
  }
}
#endif
