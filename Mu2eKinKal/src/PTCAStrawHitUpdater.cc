#include "Offline/Mu2eKinKal/inc/PTCAStrawHitUpdater.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  WireHitState PTCAStrawHitUpdater::wireHitState(ClosestApproachData const& tpdata) const {
    WireHitState whstate(WireHitState::inactive,StrawHitUpdaters::PTCA);
    if(tpdata.usable()){
      double doca = tpdata.doca();
      double absdoca = fabs(doca);
      if( absdoca < maxdoca_ && absdoca > mindoca_
          && tpdata.deltaT() < maxdt_ && tpdata.deltaT() > mindt_){  // in the sweet spot: use the DOCA to sign the ambiguity
        whstate.state_ = doca > 0.0 ? WireHitState::right : WireHitState::left;
      } else if(absdoca < mindoca_ || tpdata.deltaT() < mindt_) {
        whstate.state_ = WireHitState::null;
        whstate.nhmode_ = nhmode_;
        if(nhmode_ == WireHitState::combo){
          whstate.dvar_ =  2.0833; // (2*rstraw)^2/12   Should come from TrackerGeom TODO
        } else {
          whstate.dvar_ = dvar_;
        }
      }
    } else {
      whstate.state_ = WireHitState::unusable;
    }
    return whstate;
  }
}
