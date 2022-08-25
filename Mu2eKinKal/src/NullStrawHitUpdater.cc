#include "Offline/Mu2eKinKal/inc/NullStrawHitUpdater.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

namespace mu2e {
  using KinKal::ClosestApproachData;
  WireHitState NullStrawHitUpdater::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata) const {
    WireHitState whstate = input;
    if(input.updateable()){
      whstate.algo_ = StrawHitUpdaters::null;
      if(fabs(tpdata.doca()) < maxdoca_){
        whstate.state_ = WireHitState::null;
        whstate.nhmode_ = WireHitState::combo;
        whstate.dvar_ =  2.0833; // (2*rstraw)^2/12   Should come from TrackerGeom TODO
      } else
        whstate.state_ = WireHitState::inactive;
    }
    return whstate;
  }
}
