#ifndef TrkChargeReco_PeakFit_hh
#define TrkChargeReco_PeakFit_hh
//  Virtual class providing structure for extracting charge from ADC waveforms

#include "TrkChargeReco/inc/PeakFitParams.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  namespace TrkChargeReco {
	
    
    enum FitType {peakminusped=1,combopeakfit=2,peakfit=3};

    class PeakFit {
       
       public:
	
	// extract peak information from adc waveform data.  1 waveform generates 1 peak fit.
	// The default implementation simply sums the ADC data after subtracting pedesdal
	virtual void process(TrkTypes::ADCWaveform const& adcData, PeakFitParams & fit) const;

	// Extract charge by taking the difference of the ADC peak and pedestal
	void peakMinusPed(TrkTypes::ADCWaveform const& adcData, PeakFitParams & fit) const;

	// Initialize values to be used as parameters in fit
	void initializeFit(TrkTypes::ADCWaveform const& adcData, PeakFitParams & fit) const;

	PeakFit(const StrawElectronics& strawele, const fhicl::ParameterSet& pset);
	virtual ~PeakFit(){}


      protected:
      
	const StrawElectronics& _strawele;
	TrkChargeReco::FitType  _fittype; 

	
    };
  }
}
#endif
