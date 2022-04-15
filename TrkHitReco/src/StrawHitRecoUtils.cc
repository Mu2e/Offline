#include "Offline/TrkHitReco/inc/PeakFit.hh"
#include "Offline/TrkHitReco/inc/PeakFitRoot.hh"
#include "Offline/TrkHitReco/inc/PeakFitFunction.hh"
#include "Offline/TrkHitReco/inc/ComboPeakFitRoot.hh"

#include "Offline/DataProducts/inc/StrawEnd.hh"

#include <numeric>

#include "Offline/TrkHitReco/inc/StrawHitRecoUtils.hh"

namespace mu2e {

  void StrawHitRecoUtils::flagCrossTalk(std::unique_ptr<mu2e::StrawHitCollection> const& shCol, std::unique_ptr<mu2e::ComboHitCollection> const& chCol){
    for (size_t ilarge=0; ilarge < largeHits.size();++ilarge)
    {
      const mu2e::StrawHit& sh = (*shCol)[largeHits[ilarge]];
      for (size_t jsh : hits_by_panel[largeHitPanels[ilarge]])
      {
        if (jsh==largeHits[ilarge]) continue;
        const mu2e::StrawHit& sh2 = (*shCol)[jsh];
        if (sh2.time()-sh.time() > _ctMinT && sh2.time()-sh.time() < _ctMaxT)
        {
          if (sh.strawId().samePreamp(sh2.strawId())) (*chCol)[jsh]._flag.merge(mu2e::StrawHitFlag::elecxtalk);
          if (sh.strawId().nearestNeighbor(sh2.strawId())) (*chCol)[jsh]._flag.merge(mu2e::StrawHitFlag::strawxtalk);
        }
      }
    }
  }

  float StrawHitRecoUtils::peakMinusPedAvg(mu2e::TrkTypes::ADCWaveform const& adcData) const {
    auto wfstart = adcData.begin() + _npre;
    float pedestal = std::accumulate(adcData.begin(), wfstart, 0)*_invnpre;
    //    auto maxIter = std::max_element(wfstart,adcData.end());
    auto maxIter = wfstart;
    auto nextIter = maxIter; nextIter++;
    while(nextIter != adcData.end()){
      if (*nextIter > *maxIter)
        maxIter = nextIter;
      ++nextIter;
    }
    float peak = *maxIter;
    if(_diagLevel > 0)_maxiter->Fill(std::distance(wfstart,maxIter));
    return (peak-pedestal)*_invgainAvg;
  }

  float StrawHitRecoUtils::peakMinusPed(mu2e::StrawId id, mu2e::TrkTypes::ADCWaveform const& adcData) const {
    auto wfstart = adcData.begin() + _npre;
    float pedestal = std::accumulate(adcData.begin(), wfstart, 0)*_invnpre;
    //    auto maxIter = std::max_element(wfstart,adcData.end());
    auto maxIter = wfstart;
    while(maxIter != adcData.end() && *(maxIter+1) > *maxIter)
      ++maxIter;
    float peak = *maxIter;
    if(_diagLevel > 0)_maxiter->Fill(std::distance(wfstart,maxIter));
    return (peak-pedestal)*_invgain[id.getStraw()];
  }

  float StrawHitRecoUtils::peakMinusPedFirmware(mu2e::StrawId id, mu2e::TrkTypes::ADCValue const& pmp) const {
    return pmp*_invgain[id.getStraw()];
  }

  bool StrawHitRecoUtils::createComboHit(size_t isd, std::unique_ptr<mu2e::ComboHitCollection> const& chCol,
      std::unique_ptr<mu2e::StrawHitCollection> const& shCol, const mu2e::CaloClusterCollection* caloClusters,
      mu2e::StrawId const& sid, mu2e::TrkTypes::TDCValues const& tdc,
      mu2e::TrkTypes::TOTValues const& tot, mu2e::TrkTypes::ADCValue const& pmp, mu2e::TrkTypes::ADCWaveform const& waveform,
      mu2e::TrackerStatus const& trackerStatus, mu2e::StrawResponse const& srep, mu2e::Tracker const& tt){

    // flag digis that shouldn't be here or we don't want
    mu2e::StrawHitFlag flag;
    if (trackerStatus.noSignal(sid) || trackerStatus.suppress(sid)) {
      flag.merge(mu2e::StrawHitFlag::dead); // hits from these straws will not be used in track reconstruction
    } else if ( trackerStatus.suppress(sid)) {
      flag.merge(mu2e::StrawHitFlag::noisy); // these hits may be used in track reconstruction but not pattern recognition
    }

    // start by reconstructing the times
    mu2e::TrkTypes::TDCTimes times;
    srep.calibrateTimes(tdc,times,sid);
    // find the end with the earliest time
    mu2e::StrawEnd eend(mu2e::StrawEnd::cal);
    if(times[mu2e::StrawEnd::hv] < times[mu2e::StrawEnd::cal])
      eend = mu2e::StrawEnd(mu2e::StrawEnd::hv);
    // take the earliest of the 2 end times
    float time = times[eend.end()] - _pbtOffset;
    if (time < _minT || time > _maxT ){
      if(_filter) return false;
    } else
      flag.merge(mu2e::StrawHitFlag::timesel);

    //calorimeter filtering
    if (_usecc && caloClusters) {
      bool outsideCaloTime(true);
      for (const auto& cluster : *caloClusters)
        if (std::abs(time-cluster.time())<_clusterDt) {outsideCaloTime=false; break;}
      if (outsideCaloTime){
        if(_filter) return false;
      } else
        flag.merge(mu2e::StrawHitFlag::calosel);
    }
    // filter based on waveform shape (xtalk, undershoot, etc).  FIXME!
    //extract energy from waveform
    float energy(0.0);
    if (_fittype == mu2e::TrkHitReco::FitType::peakminuspedavg){
      float charge = peakMinusPedAvg(waveform);
      energy = srep.ionizationEnergy(charge);
    } else if (_fittype == mu2e::TrkHitReco::FitType::peakminusped){
      float charge = peakMinusPed(sid,waveform);
      energy = srep.ionizationEnergy(charge);
    } else if (_fittype == mu2e::TrkHitReco::FitType::firmwarepmp){
      float charge = peakMinusPedFirmware(sid, pmp);
      energy = srep.ionizationEnergy(charge);
    } else {
      //mu2e::TrkHitReco::PeakFitParams params;
      //_pfit->process(waveform,params);
      //energy = srep.ionizationEnergy(params._charge/srep.strawGain());
      //if (_printLevel > 1) std::cout << "Fit status = " << params._status << " NDF = " << params._ndf << " chisquared " << params._chi2
      //  << " Fit charge = " << params._charge << " Fit time = " << params._time << std::endl;
    }
    // energy selection
    if( energy > _maxE || energy < _minE ) {
      if(_filter) return false;
    } else
      flag.merge(mu2e::StrawHitFlag::energysel);
    // time-over-threshold
    mu2e::TrkTypes::TOTTimes tots{0.0,0.0};
    for(size_t iend=0;iend<2;++iend){
      tots[iend] = tot[iend]*srep.totLSB();
    }
    // choose earliest end TOT: maybe average later?
    float selected_tot = tots[eend.end()];
    // filter on specific ionization FIXME!
    // filter based on composite e/P separation FIXME!
    const mu2e::Straw& straw  = tt.getStraw( sid );
    double dw, dwerr;
    double dt = times[mu2e::StrawEnd::hv] - times[mu2e::StrawEnd::cal];
    double halfpv;
    // get distance along wire from the straw center and it's estimated error
    bool td = srep.wireDistance(straw,energy,dt, dw,dwerr,halfpv);
    float propd = straw.halfLength()+dw;
    if (eend == mu2e::StrawEnd(mu2e::StrawEnd::cal))
      propd = straw.halfLength()-dw;
    XYZVectorF pos = XYZVectorF(straw.getMidPoint()+dw*straw.getDirection());
    // create combo hit
    static const XYZVectorF _zdir(0.0,0.0,1.0);
    mu2e::ComboHit ch;
    ch._nsh = 1; // 'combo' of 1 hit
    ch._pos = pos;
    ch._wdir = straw.getDirection();
    ch._sdir = _zdir.Cross(ch._wdir);
    ch._wdist = dw;
    ch._wres = dwerr;
    ch._time = time;
    ch._edep = energy;
    ch._sid = straw.id();
    ch._dtime = srep.driftTime(straw,selected_tot,energy);
    ch._ptime = propd/(2*halfpv);
    ch._pathlength = srep.pathLength(straw,selected_tot);
    ch.addIndex(isd);
    // crude initial estimate of the transverse error
    static const float invsqrt12 = 1.0/sqrt(12.0);
    ch._tres = tt.strawOuterRadius()*invsqrt12;
    // set flags
    ch._mask = _mask;
    ch._flag = flag;
    if (td) ch._flag.merge(mu2e::StrawHitFlag::tdiv);
    ch._tend = eend;
    if(!_filter && _flagXT){
      //buffer large hit for cross-talk analysis
      size_t iplane       = straw.id().getPlane();
      size_t ipnl         = straw.id().getPanel();
      size_t global_panel = ipnl + iplane*_npanels;
      hits_by_panel[global_panel].push_back(shCol->size());
      if (energy >= _ctE) {largeHits.push_back(shCol->size()); largeHitPanels.push_back(global_panel);}
    }
    chCol->push_back(std::move(ch));
    // optionally create legacy straw hit (for diagnostics and calibration)
    if(_writesh){
      mu2e::StrawHit hit(sid,times,tots,energy);
      shCol->push_back(std::move(hit));
    }

    return true;
  }

}
