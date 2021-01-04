#ifndef TrkHitReco_PeakFit_hh
#define TrkHitReco_PeakFit_hh
//  Virtual class providing structure for extracting charge from ADC waveforms

#include "TrkHitReco/inc/PeakFitParams.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  namespace TrkHitReco {
	
    
    enum FitType {peakminuspedavg=1,peakminusped=2,combopeakfit=3,peakfit=4,firmwarepmp=5};

    class PeakFit {
       
       public:
	
	// extract peak information from adc waveform data.  1 waveform generates 1 peak fit.
	// The default implementation simply sums the ADC data after subtracting pedesdal
	virtual void process(TrkTypes::ADCWaveform const& adcData, PeakFitParams & fit) const;

	// Extract charge by taking the difference of the ADC peak and pedestal
	void peakMinusPed(TrkTypes::ADCWaveform const& adcData, PeakFitParams & fit) const;

	// Initialize values to be used as parameters in fit
	void initializeFit(TrkTypes::ADCWaveform const& adcData, PeakFitParams & fit) const;

	PeakFit(const StrawResponse& srep, const fhicl::ParameterSet& pset);
	virtual ~PeakFit(){}


      protected:
      
	const StrawResponse& _srep;
	TrkHitReco::FitType  _fittype; 

	
    };
  }
}
#endif
