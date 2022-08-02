#include "Offline/Mu2eKinKal/inc/NullStrawHitUpdater.hh"

namespace mu2e {
  using KinKal::ClosestApproachData;
  WireHitState NullStrawHitUpdater::wireHitState(ClosestApproachData const& tpdata,Straw const& straw ) const {
    // for now very simple: could check if inside straw longitudinally, etc.
    if(fabs(tpdata.doca()) > maxdoca_)
      return WireHitState(WireHitState::inactive,StrawHitUpdaters::null);
    else
      return WireHitState(WireHitState::null,StrawHitUpdaters::null);
  }
  NullHitInfo NullStrawHitUpdater::nullHitInfo(StrawResponse const& sresponse,Straw const& straw) const {
    NullHitInfo nhinfo;
    nhinfo.toff_ = dt_;
    nhinfo.dvar_ = dvar_;
    nhinfo.tvar_ = tvar_;
    nhinfo.usetime_ = usetime_;
    nhinfo.useComboDriftTime_ = true;
    return nhinfo;
  }
}
