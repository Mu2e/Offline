#include "Offline/Mu2eKinKal/inc/NullStrawHitUpdater.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

namespace mu2e {
  using KinKal::ClosestApproachData;
  WireHitState NullStrawHitUpdater::wireHitState(ClosestApproachData const& tpdata) const {
    // for now very simple: could check if inside straw longitudinally, etc.
    WireHitState whs(WireHitState::inactive,StrawHitUpdaters::null);
    if(tpdata.usable()){
      if(fabs(tpdata.doca()) < maxdoca_){
        whs.state_ = WireHitState::null;
        whs.nhmode_ = WireHitState::combo;
        whs.dvar_ =  2.0833; // (2*rstraw)^2/12   Should come from TrackerGeom TODO
      }
    } else
      whs.state_ = WireHitState::unusable;
    return whs;
  }
}
