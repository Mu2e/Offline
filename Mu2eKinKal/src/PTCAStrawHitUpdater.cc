#include "Offline/Mu2eKinKal/inc/PTCAStrawHitUpdater.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  WireHitState PTCAStrawHitUpdater::wireHitState(ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const {
    WireHitState whstate(WireHitState::inactive,StrawHitUpdaters::PTCA);
    double doca = tpdata.doca();
    double absdoca = fabs(doca);
    if( absdoca < maxdoca_ && absdoca > mindoca_
        && tpdata.deltaT() < maxdt_ && tpdata.deltaT() > mindt_){  // in the sweet spot: use the DOCA to sign the ambiguity
      whstate.state_ = doca > 0.0 ? WireHitState::right : WireHitState::left;
    } else if(absdoca < mindoca_ || tpdata.deltaT() < mindt_) {
      whstate.state_ = WireHitState::null;
    }
    // compute time and distance parameters used for null ambiguity (wire constraint)
    auto& nhinfo = whstate.nhinfo_;
    nhinfo.tmode_ = nhtmode_;
    if(nhtmode_ == NullHitInfo::usecombo){
      nhinfo.toff_ =  0.0;
      nhinfo.tvar_ = 50.0; // should come from ComboHit FIXME
      nhinfo.dvar_ = 2.1;
    } else if(nhtmode_ == NullHitInfo::usedoca){
      double vdrift = dinfo.driftVelocity_;
      nhinfo.toff_ = 0.5*mindoca_/vdrift; // this calculation is unreliable currently
      static double invthree(1.0/3.0);
      nhinfo.dvar_ = invthree*mindoca_*mindoca_;
      nhinfo.tvar_ = nhinfo.dvar_/(vdrift*vdrift);
    }
    return whstate;
  }
}
