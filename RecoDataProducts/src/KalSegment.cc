#include "RecoDataProducts/inc/KalSegment.hh"
#include "KinKal/Trajectory/CentralHelix.hh"

namespace mu2e {
  HelixVal KalSegment::helix() const {
    // CentralHelix uses the same parameter convention as BTrk.  First,
    // convert the state estimate into a helix.
    KinKal::CentralHelix chel = centralHelix();
    //convert that to the BaBar helix parameters.  We skip t0
    CLHEP::HepVector pvec(5,0);
    // crude interface to BaBar code
    pvec[KinKal::CentralHelix::d0_] = chel.d0();
    pvec[KinKal::CentralHelix::phi0_] = chel.phi0();
    pvec[KinKal::CentralHelix::omega_] = chel.omega();
    pvec[KinKal::CentralHelix::z0_] = chel.z0();
    pvec[KinKal::CentralHelix::tanDip_] = chel.tanDip();
    return HelixVal(pvec);
  }

  HelixCov KalSegment::covar() const {
    KinKal::CentralHelix chel = centralHelix();
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


  void KalSegment::mom(double flt, XYZVec& momvec) const { 
    KinKal::CentralHelix chel(_pstate,bnom());
    // momentum at the time corresponding to this flight
    auto momv = chel.momentum3(fltToTime(flt));
    // translate
    momvec = XYZVec(momv.X(),momv.Y(),momv.Z());
  }
}
