#include "Offline/Mu2eKinKal/inc/NullStrawHitUpdater.hh"

namespace mu2e {
  using KinKal::ClosestApproachData;
  WireHitState NullStrawHitUpdater::wireHitState(ClosestApproachData const& tpdata,Straw const& straw, StrawResponse const& sresponse ) const {
    // for now very simple: could check if inside straw longitudinally, etc.
    WireHitState whs;
    if(fabs(tpdata.doca()) > maxdoca_)
      whs = WireHitState(WireHitState::inactive,StrawHitUpdaters::null);
    else
      whs = WireHitState(WireHitState::null,StrawHitUpdaters::null);
    whs.nhinfo_ = NullHitInfo(0.0,tvar_,dvar_,NullHitInfo::usecombo);
    return whs;
  }
}
