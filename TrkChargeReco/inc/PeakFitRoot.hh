
#ifndef TrkChargeReco_PeakFitRoot_hh
#define TrkChargeReco_PeakFitRoot_hh
// base class for controlling peak fits
//
#include "TrkChargeReco/inc/PeakFit.hh"
#include "TrkChargeReco/inc/PeakFitFunction.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
#include <string>

class TGraphErrors;
// base class for controlling the peak fit.
namespace mu2e {

  namespace TrkChargeReco {

    //  PeakFitRoot
    //  Virtual class providing structure for extracting charge from ADC waveforms
    class PeakFitRoot : public PeakFit  {
      public:
	
	// extract peak information from adc waveform data.  1 waveform generates 1 peak fit.
	// The default implementation simply sums the ADC data after subtracting pedesdal
	virtual void process(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const;

	// Destructor
	virtual ~PeakFitRoot(){}

	// PeakFitRoot normal constructor with ConfigStruct initilization parameters
	PeakFitRoot(StrawElectronics const& strawele, FitConfig const& config, FitType const& fittype, std::string fitoptions="QNEX0S");
      	// Converts adcWaveform object to TGraphErrors object for easier manipulation in ROOT
	void adcWaveform2TGraphErrors(StrawElectronics::ADCWaveform const& adcData, TGraphErrors &fitData) const;
      protected:
	PeakFitFunction _peakfit;
	FitConfig _config;
	std::string _fitoptions;
    };
  }
}
#endif
