#include "TrkHitReco/inc/PeakFit.hh"
#include <numeric>

namespace mu2e {

   namespace TrkHitReco {


       PeakFit::PeakFit(const StrawElectronics& strawele, const fhicl::ParameterSet& pset) : 
          _strawele(strawele),
          _fittype((TrkHitReco::FitType) pset.get<TrkHitReco::FitType>("FitType",TrkHitReco::FitType::peakminusped))
       {}


       void PeakFit::process(TrkTypes::ADCWaveform const& adcData, PeakFitParams & fit) const
       {
         
         switch(_fittype)
         {
            case FitType::peakminusped :
              peakMinusPed(adcData, fit);
              break;
            default :
              initializeFit(adcData,fit);
         }
       }

       void PeakFit::peakMinusPed(TrkTypes::ADCWaveform const& adcData, PeakFitParams & fit) const
       {
          fit = PeakFitParams();
          auto maxIter = std::max_element(adcData.begin(), adcData.end());
          fit._time = std::distance(adcData.begin(), maxIter) * _strawele.adcPeriod();

          const double peak = *maxIter;
          fit._pedestal = std::accumulate(adcData.begin(), adcData.begin() + _strawele.nADCPreSamples(),0)/(double) _strawele.nADCPreSamples();

	  double charge = (peak - fit._pedestal) * _strawele.adcLSB() * _strawele.peakMinusPedestalEnergyScale();
	  fit._charge = charge;

	 //should compute width FIXME!
       }

       void PeakFit::initializeFit(const TrkTypes::ADCWaveform& adcData, PeakFitParams& fit) const
       {
          fit._earlyCharge = adcData[0]-_strawele.ADCPedestal(); // this is a crude value, should compute something FIXME!!!
          fit._pedestal = _strawele.ADCPedestal();
          fit._time = 30.0;      // this is a crude value, should compute mean or something FIXME!!!
          fit._width = 7.0;      // this is a crude value, should compute something FIXME!!!
          fit._lateShift = 50.0; // this is a crude value, should compute something FIXME!!!
          fit._lateCharge = 0.0; // 0.5*charge; this is a crude value, should compute something FIXME!!!
          fit.freeParam(PeakFitParams::charge);
       }
   }
}
