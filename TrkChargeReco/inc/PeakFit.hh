#ifndef TrkChargeReco_PeakFit_hh
#define TrkChargeReco_PeakFit_hh
// base class for controlling peak fits
//
#include "TrkChargeReco/inc/PeakFitParams.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
// base class for controlling the peak fit.
namespace mu2e {

  namespace TrkChargeReco {

    //  PeakFit
    //  Virtual class providing structure for extracting charge from ADC waveforms
    class PeakFit {
      public:
	
	// extract peak information from adc waveform data.  1 waveform generates 1 peak fit.
	// The default implementation simply sums the ADC data after subtracting pedesdal
	virtual void process(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const;

	// Destructor
	virtual ~PeakFit(){}

	// PeakFit normal constructor with ConfigStruct initilization parameters
	PeakFit(StrawElectronics const& strawele);
      private:
	StrawElectronics const& _strawele;
    };
  }
}
#endif
