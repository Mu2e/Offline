#include "Offline/Mu2eKinKal/inc/DOCAStrawHitUpdater.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  WireHitState DOCAStrawHitUpdater::wireHitState(StrawResponse const& sresponse, ClosestApproachData const& tpdata ) const {
    WireHitState whstate(WireHitState::inactive);
    double doca = tpdata.doca();
    double absdoca = fabs(doca);
    if( absdoca < maxdoca_){ // hit isn't too far from the wire
      if(absdoca > mindoca_ ){  // in the sweet spot: use the DOCA to sign the ambiguity
        whstate = doca > 0.0 ? WireHitState::right : WireHitState::left;
       // compute deltat
        if(maxdt_ > 0.0){
          auto dinfo = sresponse.driftInfoAtDistance(StrawId(),absdoca, 0.0);
          double dsign = whstate.lrSign()*tpdata.lSign(); // overall sign is the product of assigned ambiguity and doca (angular momentum) sign
          double dt = tpdata.deltaT()-dinfo.time*dsign;
          if(fabs(dt) > maxdt_) whstate = WireHitState::inactive;
        }
      } else
        whstate = WireHitState::null;
    }
    return whstate;
  }
}
