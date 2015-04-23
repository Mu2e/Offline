
#ifndef TrkChargeReco_PeakFitRoot_hh
#define TrkChargeReco_PeakFitRoot_hh
// base class for controlling peak fits
//
#include "TrkChargeReco/inc/PeakFitRootParams.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
// base class for controlling the peak fit.
namespace mu2e {

  namespace TrkChargeReco {

    //  PeakFitRoot
    //  Virtual class providing structure for extracting charge from ADC waveforms
    class PeakFitRoot : public PeakFit  {
      public:
	
	// extract peak information from adc waveform data.  1 waveform generates 1 peak fit.
	// The default implementation simply sums the ADC data after subtracting pedesdal
	virtual void process(ADCWaveform const& adcData, PeakFitRootParams & fit) const;

	// Destructor
	virtual ~PeakFitRoot(){}

	// PeakFitRoot normal constructor with ConfigStruct initilization parameters
	PeakFitRoot(StrawElectronics const& strawele, FitConfig const& config);
      private:
	PeakFitFunction _peakfit;

    };
  }
}
#endif
