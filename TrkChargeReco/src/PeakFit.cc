// MAKE SURE THAT initParams is getting passed everywhere correctly 

#include "TrkChargeReco/inc/PeakFit.hh"


//SumADC::SumADC(const ConfigStruct &initParams) : PeakFitRootBase(initParams){}
namespace mu2e {

  namespace TrkChargeReco {


    PeakFit::PeakFit(StrawElectronics const& strawele) : _strawele(strawele) {
    }

    void PeakFit::process(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const {
      double sum = 0.0;
      for (auto iadc: adcData) { sum += _strawele.adcCurrent(iadc); }
      // convert this to pC
      double charge = sum*_strawele.adcPeriod()*StrawElectronics::_pC_per_uA_ns;
      // reset
      fit = PeakFitParams();
      // set parameters
      fit._earlyCharge = adcData[0]-_strawele.ADCPedestal(); // this is a crude value, should compute something FIXME!!!
      fit._pedestal = _strawele.ADCPedestal();
      fit._time = 30.0; // this is a crude value, should compute mean or something FIXME!!!
      fit._charge = charge;
      fit._width = 7.0; // this is a crude value, should compute something FIXME!!!
      fit._lateShift = 50.0; // this is a crude value, should compute something FIXME!!!
      fit._lateCharge = 0.0; //0.5*charge; this is a crude value, should compute something FIXME!!!
      // set which parameters were free
      fit.freeParam(PeakFitParams::charge);
    }
  }
}


    
