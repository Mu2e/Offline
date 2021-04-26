#include "Mu2eKinKal/inc/KKStrawHitUpdater.hh"
//
//  Simple implementation of an ambiguity updater, provided for testing only!
//
namespace mu2e {
  void KKStrawHitUpdater::updateState(WireHitState& hitstate, ClosestApproachData const& poca) const {
    double absdoca = fabs(poca.doca()); 
    if(absdoca > mindoca_ && absdoca < maxdoca_){  // in the sweet spot: use the DOCA to sign the ambiguity
      hitstate.lrambig_ = poca.doca() > 0.0 ? WireHitState::right : WireHitState::left;
      hitstate.dimension_ = WireHitState::time;
    } else if( absdoca > maxdoca_){ // hit is too far from the wire: disable it
      hitstate.dimension_ = WireHitState::none; // disable the hit
    } else { // hit very close to the wire: ambiguity information is unusable
      hitstate.lrambig_ = WireHitState::null;
      hitstate.dimension_ = nulldim_;
    }
  }
}
