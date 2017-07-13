// MAKE SURE THAT initParams is getting passed everywhere correctly

#include "TrkChargeReco/inc/PeakFit.hh"

#include <numeric>

//SumADC::SumADC(const ConfigStruct &initParams) : PeakFitRootBase(initParams){}
namespace mu2e {

  namespace TrkChargeReco {


    PeakFit::PeakFit(StrawElectronics const& strawele, FitType const& fittype) : _strawele(strawele), _fittype(fittype) {
    }

    void PeakFit::process(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const {
      switch(_fittype) {
        case FitType::sumadc :
          sumADC(adcData, fit);
          break;
        case FitType::peakminusped :
          peakMinusPed(adcData, fit);
          break;
        default :
          // set parameters for fitting
          initializeFit(adcData,fit);
      }
    }

    void PeakFit::sumADC(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const{
      // reset
      fit = PeakFitParams();
      double sum = 0.0;
      for (auto iadc : adcData) { sum += _strawele.adcCurrent(iadc); }
      // convert charge to pC
      double charge = sum*_strawele.adcPeriod()*StrawElectronics::_pC_per_uA_ns;
      fit._charge = charge;
    }

    void PeakFit::peakMinusPed(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const
    {
      // reset
      fit = PeakFitParams();
      auto maxIter = std::max_element(adcData.begin(), adcData.end());
      fit._time = std::distance(adcData.begin(), maxIter) * _strawele.adcPeriod();

      const double peak = *maxIter;
      double pedestal = std::accumulate(adcData.begin(), adcData.begin() + _strawele.nADCPreSamples(),0)/(double) _strawele.nADCPreSamples();

      // convert charge to pC
     double charge = (peak - pedestal) * _strawele.adcLSB() / _strawele.normalization(StrawElectronics::adc) / exp(-1.0) * _strawele.peakMinusPedestalEnergyScale();
     fit._charge = charge;
    }

    void PeakFit::initializeFit(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const{
      // set parameters
      fit._earlyCharge = adcData[0]-_strawele.ADCPedestal(); // this is a crude value, should compute something FIXME!!!
      fit._pedestal = _strawele.ADCPedestal();
      fit._time = 30.0; // this is a crude value, should compute mean or something FIXME!!!
      fit._width = 7.0; // this is a crude value, should compute something FIXME!!!
      fit._lateShift = 50.0; // this is a crude value, should compute something FIXME!!!
      fit._lateCharge = 0.0; //0.5*charge; this is a crude value, should compute something FIXME!!!
      // set which parameters were free
      fit.freeParam(PeakFitParams::charge);
    }
  }
}
