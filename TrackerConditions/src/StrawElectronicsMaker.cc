
#include "Offline/TrackerConditions/inc/StrawElectronicsMaker.hh"
#include "Offline/GeneralUtilities/inc/DigitalFiltering.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "cetlib_except/exception.h"
#include "TMath.h"
#include <cmath>
#include <complex>
#include <memory>

using namespace std;

namespace mu2e {
  using namespace TrkTypes;

  namespace {
    // a per-channel fcl list is either empty (no override) or covers every channel
    void checkChannelListLength(char const* name, size_t size, size_t nchannels) {
      if (size != 0 && size != nchannels)
        throw cet::exception("BADCONFIG")
          << "StrawElectronics fcl parameter " << name << " has " << size
          << " entries; it must be empty or have " << nchannels << "\n";
    }
  }

  StrawElectronics::ptr_t StrawElectronicsMaker::fromFcl(EventTiming::cptr_t eventTiming) {

    unsigned maxTDC = (0x1<<_config.numTDCbits())-1;

    // the wire distance points are interpolated between neighbours, so there must be at least 2
    auto const wireDistances = _config.wireDistances();
    auto const currentMeans = _config.currentMeans();
    auto const currentNormalizations = _config.currentNormalizations();
    auto const currentSigmas = _config.currentSigmas();
    auto const currentT0s = _config.currentT0s();
    if (wireDistances.size() < 2
        || currentMeans.size() != wireDistances.size()
        || currentNormalizations.size() != wireDistances.size()
        || currentSigmas.size() != wireDistances.size()
        || currentT0s.size() != wireDistances.size()) {
      throw cet::exception("BADCONFIG")
        << "StrawElectronics fcl parameters wireDistances, currentMeans, currentNormalizations,"
        << " currentSigmas and currentT0s must have the same length, at least 2; lengths are "
        << wireDistances.size() << " " << currentMeans.size() << " " << currentNormalizations.size()
        << " " << currentSigmas.size() << " " << currentT0s.size() << "\n";
    }

    // creat this at the beginning since it must be used,
    // partially constructed, to complete the construction
    auto ptr = std::make_shared<StrawElectronics>(
        _config.overrideDbTimeOffsets(), _config.deadTimeAnalog(),
        _config.deadTimeDigital(), _config.saturationVoltage(), _config.strawNoise(),
        _config.ADCLSB(), _config.maxADC(), _config.nADCPackets(), _config.nADCPresamples(),
        _config.ADCPeriod(), _config.ADCOffset(),
        _config.maxThreshTimeSeparation(),
        _config.TDCLSB(), maxTDC, _config.TOTLSB(),
        _config.maxTOT(), _config.TDCResolution(),
        _config.electronicsTimeDelay(), _config.eventWindowMarkerROCJitter(),
        _config.digitizationStart(),
        _config.digitizationEnd(),
        _config.responseBins(),
        _config.sampleRate(), _config.saturationSampleFactor(),
        _config.preampPoles(), _config.preampZeros(),
        _config.ADCPoles(), _config.ADCZeros(),
        _config.preampToAdc1Poles(), _config.preampToAdc1Zeros(),
        _config.preampToAdc2Poles(), _config.preampToAdc2Zeros(),
        wireDistances, currentMeans,
        currentNormalizations, currentSigmas,
        currentT0s, _config.reflectionTimeShift(),
        _config.reflectionVelocity(), _config.reflectionALength(),
        _config.reflectionFrac(), _config.triggerHysteresis(),
        _config.clusterLookbackTime());

    auto const timeOffsetPanelFcl = _config.timeOffsetPanel();
    auto const timeOffsetStrawHVFcl = _config.timeOffsetStrawHV();
    auto const timeOffsetStrawCalFcl = _config.timeOffsetStrawCal();
    checkChannelListLength("timeOffsetPanel", timeOffsetPanelFcl.size(), StrawId::_nupanels);
    checkChannelListLength("timeOffsetStrawHV", timeOffsetStrawHVFcl.size(), StrawId::_nustraws);
    if (timeOffsetStrawCalFcl.size() != timeOffsetStrawHVFcl.size())
      throw cet::exception("BADCONFIG")
        << "StrawElectronics fcl parameters timeOffsetStrawHV and timeOffsetStrawCal must have the same length, not "
        << timeOffsetStrawHVFcl.size() << " and " << timeOffsetStrawCalFcl.size() << "\n";
    std::array<double, StrawId::_nupanels> timeOffsetPanel;
    std::array<double, StrawId::_nustraws> timeOffsetStrawHV, timeOffsetStrawCal;
    if (timeOffsetPanelFcl.size() > 0){
      for (size_t i=0;i<timeOffsetPanel.size();i++)
        timeOffsetPanel[i] = timeOffsetPanelFcl[i];
    }else{
      timeOffsetPanel.fill(0);
    }
    if (timeOffsetStrawHVFcl.size() > 0){
      for(size_t i=0;i<timeOffsetStrawHV.size();i++){
        timeOffsetStrawHV[i] = timeOffsetStrawHVFcl[i];
        timeOffsetStrawCal[i] = timeOffsetStrawCalFcl[i];
      }
    }else{
      timeOffsetStrawHV.fill(0);
      timeOffsetStrawCal.fill(0);
    }
    ptr->setOffsets( timeOffsetPanel,
        timeOffsetStrawHV,
        timeOffsetStrawCal );

    auto const thresholddVdI = _config.thresholddVdI();
    auto const adcdVdI = _config.adcdVdI();
    checkChannelListLength("thresholddVdI", thresholddVdI.size(), StrawId::_nustraws);
    checkChannelListLength("adcdVdI", adcdVdI.size(), StrawId::_nustraws);
    std::array<std::array<double, StrawId::_nustraws>,StrawElectronics::npaths> dVdI;
    if(thresholddVdI.size()>0) {
      for (size_t i=0;i<dVdI[StrawElectronics::thresh].size();i++)
        dVdI[StrawElectronics::thresh][i] = thresholddVdI[i];
    } else {
      dVdI[StrawElectronics::thresh].fill(_config.defaultThresholddVdI());
    }
    if(adcdVdI.size()>0) {
      for (size_t i=0;i<dVdI[StrawElectronics::adc].size();i++)
        dVdI[StrawElectronics::adc][i] = adcdVdI[i];
    } else {
      dVdI[StrawElectronics::adc].fill(_config.defaultAdcdVdI());
    }
    ptr->setdVdI(dVdI[StrawElectronics::thresh],StrawElectronics::thresh);
    ptr->setdVdI(dVdI[StrawElectronics::adc],StrawElectronics::adc);

    std::array<double,StrawElectronics::npaths> analognoise;
    analognoise[StrawElectronics::thresh] = _config.thresholdAnalogNoise();
    analognoise[StrawElectronics::adc]    = _config.adcAnalogNoise();
    ptr->setAnalogNoise(analognoise);

    auto const discriminatorThreshold = _config.discriminatorThreshold();
    checkChannelListLength("discriminatorThreshold", discriminatorThreshold.size(), StrawId::_nustrawends);
    std::array<double, StrawId::_nustrawends> vthresh;
    if(discriminatorThreshold.size()>0) {
      for (size_t i=0;i<vthresh.size();i++)
        vthresh[i] = discriminatorThreshold[i];
    } else {
      vthresh.fill(_config.defaultDiscriminatorThreshold());
    }
    ptr->setvthresh(vthresh);

    // here we start using the partially constructed StrawElectronics *ptr
    ptr->setDigitizationStartTDC( ptr->tdcResponse( _config.digitizationStart() - _config.electronicsTimeDelay() - eventTiming->timeFromProtonsToDRMarker() ));
    ptr->setDigitizationEndTDC( ptr->tdcResponse( _config.digitizationEnd() - _config.electronicsTimeDelay() - eventTiming->timeFromProtonsToDRMarker()  ));
    ptr->setTimeFromProtonsToDRMarker(eventTiming->timeFromProtonsToDRMarker());

    std::array<uint16_t, StrawId::_nustraws> ADCped;
    for (size_t i=0;i<ADCped.size();i++){
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

    auto integral_normalization = log(responseBins/2/sampleRate + currentT0s[0]) - log(currentT0s[0]); // integral of 1/(t+t0) for 0 cm

    std::vector<StrawElectronics::WireDistancePoint> wPoints;
    for (size_t ai=0;ai<wireDistances.size();ai++){
      wPoints.emplace_back( wireDistances[ai],
          currentMeans[ai],  currentNormalizations[ai],
          currentSigmas[ai], currentT0s[ai]);
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
        wPoints[ai]._currentPulse[i] /= integral_normalization; // normalization for 1/(t+t0)
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
      double preampToAdc1Max = 0;
      double linmax_thresh = 0;
      double linmax_adc = 0;
      for (int i=0;i<responseBins;i++){
        if (wPoints[ai]._preampResponse[i] > linmax_thresh){
          linmax_thresh = wPoints[ai]._preampResponse[i];
        }
        if (wPoints[ai]._adcResponse[i] > linmax_adc){
          linmax_adc = wPoints[ai]._adcResponse[i];
        }

        if (wPoints[ai]._preampToAdc1Response[i] > preampToAdc1Max)
          preampToAdc1Max = wPoints[ai]._preampToAdc1Response[i];
      }

      // normalize preampToAdc1Response to match preampResponse
      for (int i=0;i<responseBins;i++){
        wPoints[ai]._preampToAdc1Response[i] *= linmax_thresh/preampToAdc1Max;
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

    double linmax_adc = 0;
    for (int i=0;i<responseBins;i++){
      if (wPoints[0]._adcResponse[i] > linmax_adc){
        linmax_adc = wPoints[0]._adcResponse[i];
      }
    }

    for (int i=0;i<responseBins;i++){
      preampToAdc2Response[i] *= linmax_adc/preampToAdc2Max;
      preampToAdc2test[i] *= linmax_adc/preampToAdc2Max;
    }

    ptr->setPreampToAdc2Response(preampToAdc2Response);

    GeomHandle<mu2e::Tracker> th;
    const Tracker* tracker = th.get();

    for (size_t ipoint=0;ipoint<wPoints.size();ipoint++){
      for (size_t ipath=0;ipath<StrawElectronics::npaths;ipath++){
        for (size_t is=0;is<StrawId::_nstraws;is++){
          StrawId sid(0,0,is);
          Straw const& straw = tracker->getStraw(sid);
          double tmax = 0;
          double linmax = -9e9;
          int start = -5*sampleRate;
          int end = 20*sampleRate;
          for (int k=start;k<end;k++){
            double time = k/(double)sampleRate;
            double resp = ptr->linearResponse(straw,static_cast<StrawElectronics::Path>(ipath),time,1e-9,wPoints[ipoint]._distance,false);
            resp /= dVdI[ipath][sid.uniqueStraw()]; // linmax should be normalized to dVdI = 1
            if (resp > linmax){
              linmax = resp;
              tmax = time;
            }
          }
          wPoints[ipoint]._tmax[ipath][is] = tmax;
          wPoints[ipoint]._linmax[ipath][is] = linmax*1e9;
        }
      }
    }

    ptr->setwPoints(wPoints);
    return ptr;

  } // end fromFcl

  StrawElectronics::ptr_t StrawElectronicsMaker::fromDb(
      TrkDelayPanel::cptr_t tdp,
      TrkDelayRStraw::cptr_t tdrs,
      TrkPreampStraw::cptr_t tps,
      EventTiming::cptr_t eventTiming ) {
    // initially fill from fcl to get all the constants
    auto ptr = fromFcl(eventTiming);

    // now overwrite with db values

    std::array<double, StrawId::_nustrawends> vthresh;
    for(size_t i=0; i<StrawId::_nustraws; i++) {
      vthresh[2*i+StrawEnd::cal] = tps->rowAt(i).thresholdCal();
      vthresh[2*i+StrawEnd::hv] = tps->rowAt(i).thresholdHv();
    }

    ptr->setvthresh(vthresh);

    std::array<uint16_t, StrawId::_nustraws> ADCped;
    for (size_t i=0;i<ADCped.size();i++){
      double avgThresh = ( vthresh[i*2+0] +  vthresh[i*2+1])/2.;
      ADCped[i] = (ADCValue) ((_config.maxADC()+1)/2. - avgThresh/_config.ADCLSB() * 20);
    }
    ptr->setADCPed(ADCped);

    std::array<double, StrawId::_nupanels> timeOffsetPanel;
    std::array<double, StrawId::_nustraws> timeOffsetStrawHV;
    std::array<double, StrawId::_nustraws> timeOffsetStrawCal;

    for(size_t i=0; i<StrawId::_nupanels; i++) {
      timeOffsetPanel[i] = tdp->rowAt(i).delay();
    }
    for(size_t i=0; i<StrawId::_nustraws; i++) {
      size_t istraw = i % StrawId::_nstraws;
      timeOffsetStrawHV[i] = tdrs->rowAt(istraw).delayHv() + tps->rowAt(i).delayHv();
      timeOffsetStrawCal[i] = tdrs->rowAt(istraw).delayCal() + tps->rowAt(i).delayCal();
    }

    ptr->setOffsets( timeOffsetPanel,
        timeOffsetStrawHV,
        timeOffsetStrawCal );

    std::array<double, StrawId::_nustraws> adcdVdI;
    for (size_t i=0;i<StrawId::_nustraws;i++) {
      adcdVdI[i] = tps->rowAt(i).gain();
    }
    ptr->setdVdI(adcdVdI,StrawElectronics::adc);

    if (_config.verbose())
      ptr->print(std::cout);

    return ptr;

  } // end fromDb

}
