#ifndef TrkChargeRecoOld_PeakFit_hh
#define TrkChargeRecoOld_PeakFit_hh
// base class for controlling peak fits
//
#include "TrkChargeRecoOld/inc/PeakFitParams.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
// base class for controlling the peak fit.
namespace mu2e {

  namespace TrkChargeRecoOld {
	
	enum FitType {sumadc=0,peakminusped=1,combopeakfit=2,peakfit=3};

    //  PeakFit
    //  Virtual class providing structure for extracting charge from ADC waveforms
    class PeakFit {
      public:
	
	// extract peak information from adc waveform data.  1 waveform generates 1 peak fit.
	// The default implementation simply sums the ADC data after subtracting pedesdal
	virtual void process(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const;

	// Extract charge by summing the ADC data after subtracting pedestal
	void sumADC(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const;
	
	// Extract charge by taking the difference of the ADC peak and pedestal
	void peakMinusPed(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const;

	// Initialize values to be used as parameters in fit
	void initializeFit(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const;

	// Destructor
	virtual ~PeakFit(){}

	// PeakFit normal constructor with ConfigStruct initilization parameters
	PeakFit(StrawElectronics const& strawele, FitType const& fittype);

      protected:
	StrawElectronics const& _strawele;

	FitType const& _fittype; 

	
    };
  }
}
#endif
