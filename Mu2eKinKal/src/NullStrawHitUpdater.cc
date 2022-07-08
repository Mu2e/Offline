#include "Offline/Mu2eKinKal/inc/NullStrawHitUpdater.hh"
namespace mu2e {
  using KinKal::ClosestApproachData;
  WireHitState NullStrawHitUpdater::wireHitState(ClosestApproachData const& tpdata ) const {
    if(fabs(tpdata.doca()) > maxdoca_)
      return WireHitState(WireHitState::inactive);
    else
      return WireHitState(WireHitState::null);
  }
}
