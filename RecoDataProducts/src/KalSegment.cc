#include "Offline/RecoDataProducts/inc/KalSegment.hh"

namespace mu2e {

  KinKal::TimeRange KalSegment::timeRange() const {
    // protect against unphysical times
    if(_tmin <= _tmax){
      return KinKal::TimeRange(_tmin,_tmax);
    } else {
      std::cout << "KalSegment: Invalid time range, tmin " << _tmin << " tmax " << _tmax << std::endl;
      return KinKal::TimeRange(_pstate.time(),_pstate.time());
    }
  }

  double KalSegment::t0Val(TrkFitFlag const& flag) const {
    // convert to the appropriate parameterization.
    if(flag.hasAllProperties(TrkFitFlag::KKLoopHelix)) {
      auto lh = loopHelix();
      return lh.paramVal(KinKal::LoopHelix::t0Index());
    } else if(flag.hasAllProperties(TrkFitFlag::KKCentralHelix)) {
      auto ch = loopHelix();
      return ch.paramVal(KinKal::CentralHelix::t0Index());
    } else if (flag.hasAllProperties(TrkFitFlag::KKLine)) {
      auto kl = kinematicLine();
      return kl.paramVal(KinKal::KinematicLine::t0Index());
    }
    //      throw cet::exception("RECO")<<"mu2e::KalSegment: no trajectory specified in flag" << std::endl;
    // for now, revert to the legacy implementation.  Once BTrk is fully removed this should be removed
    auto vel = _pstate.velocity();
    return _pstate.time() - _pstate.position3().Z()/vel.Z();
  }

  // deprecated legacy functions, these should be removed FIXME

  HelixVal KalSegment::helix() const {
    // CentralHelix uses the same parameter convention as BTrk.  First,
    // convert the state estimate into a helix.
    auto chel = centralHelix();
    //convert that to the BaBar helix parameters.  skip t0 since that isn't part of the geometric helix
    CLHEP::HepVector pvec(5,0);
    pvec[KinKal::CentralHelix::d0_] = chel.d0();
    pvec[KinKal::CentralHelix::phi0_] = chel.phi0();
    pvec[KinKal::CentralHelix::omega_] = chel.omega();
    pvec[KinKal::CentralHelix::tanDip_] = chel.tanDip();
    pvec[KinKal::CentralHelix::z0_] = chel.z0();
    return HelixVal(pvec);
  }

  HelixCov KalSegment::covar() const {
    auto chel = centralHelix();
    auto const& kkcov = chel.params().covariance();
    CLHEP::HepSymMatrix cov(5,1);
    for(size_t ipar=0; ipar <5; ipar++){
      for(size_t jpar=0; jpar <=ipar; jpar++){
        // stupid fotran-like interface
        cov.fast(ipar+1,jpar+1) = kkcov(ipar,jpar);
      }
    }
    return HelixCov(cov);
  }

  void KalSegment::mom(double flt, XYZVectorF& momvec) const {
    auto chel = centralHelix();
    // momentum at the time corresponding to this flight
    auto momv = chel.momentum3(fltToTime(flt));
    // translate
    momvec = XYZVectorF(momv.X(),momv.Y(),momv.Z());
  }

  double KalSegment::fltToTime(double flt) const {
    auto vel = _pstate.velocity();
    return t0Val() + flt/vel.R();
  }

  double KalSegment::timeToFlt(double time) const {
    auto vel = _pstate.velocity();
    return (time - t0Val())*vel.R();
  }
}
