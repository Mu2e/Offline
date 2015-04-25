// fit waveform using root TF1
#include "TrkChargeReco/inc/PeakFitRoot.hh"
#include "TF1.h"
#include "TGraphErrors.h"

namespace mu2e {
  namespace TrkChargeReco {

  PeakFitRoot::PeakFitRoot(StrawElectronics const& strawele, FitConfig const& config) : PeakFit(strawele),
  _peakfit(strawele,config) {}

  void PeakFitRoot::process(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const {
    // find initial values for the fit
    PeakFit::process(adcData,fit);
    // set the initial values based on this 'fit'
    Double_t parray[PeakFitParams::nParams];
    fit.fillArray(parray);
    _peakfit.fitModelTF1()->SetParameters(parray);
    // cobnvert waveform to a TGraph
    TGraphErrors fitData;
    adcWaveform2TGraphErrors(adcData,fitData);
    // invoke the fit
    fitData.Fit(_peakfit.fitModelTF1(),"QN"); // fit options should be configurable, FIXME!!
    // copy fit result back to the the PeakFitParams
    fit = PeakFitParams(_peakfit.fitModelTF1()->GetParameters());
  }

  void PeakFitRoot::adcWaveform2TGraphErrors(StrawElectronics::ADCWaveform const& adcData, TGraphErrors &fitData) const{
      Double_t adcDataTemp[adcData.size()];
      Double_t measurementTimes[adcData.size()];
      Double_t measurementTimesErrors[adcData.size()];
      Double_t adcDataErrors[adcData.size()];

      for (size_t i = 0; i < adcData.size(); ++i)
      {
	adcDataTemp[i] = (Double_t) adcData[i];
	measurementTimes[i] = (Double_t) i * _strawele.adcPeriod(); // should deal with global time offset FIXME!!
	measurementTimesErrors[i] = _strawele.adcPeriod();
	adcDataErrors[i] = 1.0*_strawele.analogNoise(StrawElectronics::adc)/_strawele.adcLSB(); // should be able to scale the error FIXME!!
      }
      fitData = TGraphErrors(adcData.size(),measurementTimes,adcDataTemp,measurementTimesErrors,adcDataErrors);
    }

  } // TrkChargeReco namespace

}// mu2e namespace


