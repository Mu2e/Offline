
#include "TrackerConditions/inc/StrawElectronicsMaker.hh"
#include "GeneralUtilities/inc/DigitalFiltering.hh"
#include "cetlib_except/exception.h"
#include "TMath.h"
#include <math.h>
#include <complex>
#include <memory>

using namespace std;

namespace mu2e {
  using namespace TrkTypes;

  StrawElectronics::ptr_t StrawElectronicsMaker::fromFcl() {

    // creat this at the beginning since it must be used,
    // partially constructed, to complete the construction
    auto ptr = std::make_shared<StrawElectronics>(_config.deadTimeAnalog(), 
       _config.deadTimeDigital(), _config.saturationVoltage(), _config.strawNoise(), 
       _config.ADCLSB(), _config.maxADC(), _config.nADCPresamples(),
       _config.ADCPeriod(), _config.ADCOffset(), 
       _config.maxThreshTimeSeparation(), _config.tCoince(),
       _config.TDCLSB(), _config.maxTDC(), _config.TOTLSB(), 
       _config.maxTOT(), _config.TDCResolution(), 
       _config.electronicsTimeDelay(), _config.eventWindowMarkerROCJitter(), 
       _config.flashStart(), _config.flashEnd(), _config.flashClockSpeed(),
       _config.responseBins(), 
       _config.sampleRate(), _config.saturationSampleFactor(), 
       _config.preampPoles(), _config.preampZeros(), 
       _config.ADCPoles(), _config.ADCZeros(), 
       _config.preampToAdc1Poles(), _config.preampToAdc1Zeros(), 
       _config.preampToAdc2Poles(), _config.preampToAdc2Zeros(), 
       _config.wireDistances(), _config.currentMeans(), 
       _config.currentNormalizations(), _config.currentSigmas(), 
       _config.currentT0s(), _config.clusterLookbackTime(), 
       _config.timeOffsetPanel(), _config.timeOffsetStrawHV(), 
       _config.timeOffsetStrawCal() );

    std::array<vector<double>,StrawElectronics::npaths> dVdI;
    if(_config.thresholddVdI().size()>0) {
      dVdI[StrawElectronics::thresh] = _config.thresholddVdI();
    } else {
      dVdI[StrawElectronics::thresh] = vector<double>(96,_config.defaultThresholddVdI());
    }
    if(_config.adcdVdI().size()>0) {
      dVdI[StrawElectronics::adc] = _config.adcdVdI();
    } else {
      dVdI[StrawElectronics::adc] = vector<double>(96,_config.defaultAdcdVdI());
    }
    ptr->setdVdI(dVdI);

    std::array<double,StrawElectronics::npaths> analognoise;
    analognoise[StrawElectronics::thresh] = _config.thresholdAnalogNoise();
    analognoise[StrawElectronics::adc]    = _config.adcAnalogNoise();
    ptr->setAnalogNoise(analognoise);
    
    vector<double> vthresh;
    if(_config.discriminatorThreshold().size()>0) {
      vthresh = _config.discriminatorThreshold();
    } else {
      vthresh = vector<double>(192,_config.defaultDiscriminatorThreshold());
    }
    ptr->setvthresh(vthresh);

    // here we start using the partially constructed StrawElectronics *ptr
    ptr->setFlashTDC( ptr->tdcResponse( _config.flashStart() ),
		      ptr->tdcResponse( _config.flashEnd()   ) );

    std::vector<uint16_t> ADCped(96,0);
    for (int i=0;i<96;i++){
      double avgThresh = ( vthresh[i*2+0] +  vthresh[i*2+1])/2.;
      ADCped[i] = (ADCValue) ((_config.maxADC()+1)/2. - avgThresh/_config.ADCLSB() * 20);
    }
    ptr->setADCPed(ADCped);

    std::array<double,StrawElectronics::npaths> ttrunc;
    auto responseBins = _config.responseBins();
    auto sampleRate = _config.sampleRate();
    ttrunc[StrawElectronics::thresh] = (responseBins/2)/_config.sampleRate();
    ttrunc[StrawElectronics::adc] = (responseBins/2)/_config.sampleRate();
    ptr->setttrunc(ttrunc);

    // precompute ion drift current pulses
    auto currentImpulse = std::vector<double>(responseBins,0);
    currentImpulse[responseBins/2] = 1;
    auto preampToAdc2Response = std::vector<double>(responseBins,0);

    ptr->calculateResponse(_config.preampToAdc2Poles(),_config.preampToAdc2Zeros(),
			   currentImpulse,preampToAdc2Response);

    ptr->setCurrentImpulse(currentImpulse);

    const double pC_per_uA_ns{1000}; // unit conversion from pC/ns to microAmp
    auto thresh = StrawElectronics::thresh;
    auto adc    = StrawElectronics::adc;

    std::vector<StrawElectronics::WireDistancePoint> wPoints;
    for (size_t ai=0;ai<_config.wireDistances().size();ai++){
      wPoints.emplace_back( _config.wireDistances()[ai],
	       _config.currentMeans()[ai],  _config.currentNormalizations()[ai],
               _config.currentSigmas()[ai], _config.currentT0s()[ai]);
      wPoints[ai]._currentPulse = std::vector<double>(responseBins,0);
      double integral = 0;
      for (int i=0;i<responseBins;i++){
	double t_gaus = (i-responseBins/2)/sampleRate;
	double val_gaus = 1/sqrt(TMath::TwoPi()*wPoints[ai]._sigma*wPoints[ai]._sigma)*exp(-((t_gaus-wPoints[ai]._mean)*(t_gaus-wPoints[ai]._mean))/(2*wPoints[ai]._sigma*wPoints[ai]._sigma));
	for (int j=0;j<responseBins-i;j++){
	  double t_tail = j/sampleRate;
	  double val = val_gaus / (t_tail + wPoints[ai]._t0);
	  wPoints[ai]._currentPulse[i+j] += val;
	  integral += val;
	}
      }
      for (int i=0;i<responseBins;i++){
	// correct for sampleRate so that calculateResponse peak is independent of it
	// this combined with pC_per_uA_ns is the unit transform from pC to uA
	wPoints[ai]._currentPulse[i] /= sampleRate * pC_per_uA_ns;
	wPoints[ai]._currentPulse[i] /= 4.615; // normalization for 1/(t+t0)
	wPoints[ai]._currentPulse[i] *= wPoints[ai]._normalization;
      }
      
      // calculate parameters for transfer function
      wPoints[ai]._preampResponse = std::vector<double>(responseBins,0);
      wPoints[ai]._adcResponse = std::vector<double>(responseBins,0);
      wPoints[ai]._preampToAdc1Response = std::vector<double>(responseBins,0);
      ptr->calculateResponse(_config.preampPoles(),_config.preampZeros(),
			wPoints[ai]._currentPulse,wPoints[ai]._preampResponse);
      ptr->calculateResponse(_config.ADCPoles(),_config.ADCZeros(),
			wPoints[ai]._currentPulse,wPoints[ai]._adcResponse);
      ptr->calculateResponse(_config.preampToAdc1Poles(),_config.preampToAdc1Zeros(),
			wPoints[ai]._currentPulse,wPoints[ai]._preampToAdc1Response);
      
      // now set other parameters
      wPoints[ai]._tmax[thresh] = 0;
      wPoints[ai]._linmax[thresh] = 0;
      wPoints[ai]._tmax[adc] = 0;
      wPoints[ai]._linmax[adc] = 0;
      double preampToAdc1Max = 0;
      for (int i=0;i<responseBins;i++){
	if (wPoints[ai]._preampResponse[i] > wPoints[ai]._linmax[thresh]){
	  wPoints[ai]._linmax[thresh] = wPoints[ai]._preampResponse[i];
	  wPoints[ai]._tmax[thresh] = (-responseBins/2 + i)/sampleRate;
	}
	if (wPoints[ai]._adcResponse[i] > wPoints[ai]._linmax[adc]){
	  wPoints[ai]._linmax[adc] = wPoints[ai]._adcResponse[i];
	  wPoints[ai]._tmax[adc] = (-responseBins/2 + i)/sampleRate;
	}
	if (wPoints[ai]._preampToAdc1Response[i] > preampToAdc1Max)
	  preampToAdc1Max = wPoints[ai]._preampToAdc1Response[i];
      }
      
      // normalize preampToAdc1Response to match preampResponse
      for (int i=0;i<responseBins;i++){
	wPoints[ai]._preampToAdc1Response[i] *= 
	  wPoints[ai]._linmax[thresh]/preampToAdc1Max;
      }
    } // loop over wireDistances

    // normalize preampToAdc2Response to match adcResponse
    std::vector<double> preampToAdc2test(responseBins,0);
    ptr->calculateResponse(_config.preampToAdc2Poles(),_config.preampToAdc2Zeros(),
		      wPoints[0]._preampToAdc1Response,preampToAdc2test);

    ptr->setwPoints(wPoints);

    double preampToAdc2Max = 0;
    for (int i=0;i<responseBins;i++){
      if (preampToAdc2test[i] > preampToAdc2Max)
	preampToAdc2Max = preampToAdc2test[i];
    }

    for (int i=0;i<responseBins;i++){
      preampToAdc2Response[i] *= wPoints[0]._linmax[adc]/preampToAdc2Max;
      preampToAdc2test[i] *= wPoints[0]._linmax[adc]/preampToAdc2Max;
    }
    
    ptr->setPreampToAdc2Response(preampToAdc2Response);

    return ptr;

  } // end fromFcl

  StrawElectronics::ptr_t StrawElectronicsMaker::fromDb(
				   TrkDelayPanel::cptr_t tdp,
				   TrkPreampRStraw::cptr_t tprs,
				   TrkPreampStraw::cptr_t tps,
				   TrkThresholdRStraw::cptr_t ttrs ) {
    // initially fill from fcl to get all the constants
    auto ptr = fromFcl();

    // now overwrite with db values

    vector<double> vthresh(2*StrawId::_nstraws);
    for(size_t i=0; i<StrawId::_nstraws; i++) {
      vthresh[2*i+StrawEnd::cal] = ttrs->rowAt(i).thresholdCal();
      vthresh[2*i+StrawEnd::hv] = ttrs->rowAt(i).thresholdHv();
    }

    ptr->setvthresh(vthresh);


    std::vector<double> timeOffsetPanel(StrawId::_nupanels);
    std::vector<double> timeOffsetStrawHV(StrawId::_nstraws);
    std::vector<double> timeOffsetStrawCal(StrawId::_nstraws);
    for(size_t i=0; i<StrawId::_nupanels; i++) {
      timeOffsetPanel[i] = tdp->rowAt(i).delay();
    }
    for(size_t i=0; i<StrawId::_nstraws; i++) {
      timeOffsetStrawHV[i] = tprs->rowAt(i).delayHv();
      timeOffsetStrawCal[i] = tprs->rowAt(i).delayCal();
    }

    ptr->setOffsets( timeOffsetPanel,
		     timeOffsetStrawHV,
		     timeOffsetStrawCal );

    return ptr;

  } // end fromDb

}
