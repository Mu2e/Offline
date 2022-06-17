#include "Offline/Mu2eKinKal/inc/DOCAStrawHitUpdater.hh"
#include <cmath>

namespace mu2e {
  KinKal::WireHitState DOCAStrawHitUpdater::wireHitState(double doca) const {
    using KinKal::WireHitState;
    WireHitState whstate(WireHitState::inactive);
    double absdoca = fabs(doca);
    if( absdoca < maxdoca_){ // hit isn't too far from the wire
      if(absdoca > mindoca_ && absdoca < maxddoca_){  // in the sweet spot: use the DOCA to sign the ambiguity
        whstate = doca > 0.0 ? WireHitState::right : WireHitState::left;
      } else { // hit too close to the wire to resolve ambiguity, or with a suspiciously large drift: just use the raw wire position and time to constrain the track
        whstate = WireHitState::null;
      }
    }
    return whstate;
  }
}
