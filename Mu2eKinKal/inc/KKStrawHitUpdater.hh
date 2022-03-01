#ifndef Mu2eKinKal_KKStrawHitUpdater_hh
#define Mu2eKinKal_KKStrawHitUpdater_hh

#include "KinKal/Detector/WireHitStructs.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"

namespace mu2e {
  // interface for updating individual straw hits
  template <class KTRAJ> class KKStrawHit;
  class KKStrawHitUpdater {
    public:
      KKStrawHitUpdater() : maxdoca_(-1.0), minprob_(-1.0), mindoca_(0), maxddoca_(-1) {}
      KKStrawHitUpdater(double maxdoca, double minprob, double minddoca, double maxddoca) :
        maxdoca_(maxdoca), minprob_(minprob),
        mindoca_(minddoca), maxddoca_(maxddoca) {}
      template <class KTRAJ> void update(KKStrawHit<KTRAJ>& swh) const;
    private:
      double maxdoca_; // maximum DOCA to still use a hit
      double minprob_; // minimum chisqquared probability to use a hit
      double mindoca_; // minimum DOCA to use drift information
      double maxddoca_; // maximum DOCA to use drift information
  };

  template <class KTRAJ> void KKStrawHitUpdater::update(KKStrawHit<KTRAJ>& kksh) const {
    using KinKal::ClosestApproachData;
    using KinKal::WireHitState;
    kksh.mindoca_ = std::min(mindoca_,kksh.strawRadius());
    WireHitState::State state;
    if(kksh.closestApproach().usable()){
      auto doca = kksh.closestApproach().doca();
      auto absdoca = fabs(doca);
      auto chisq = kksh.chisq();
      if( absdoca > maxdoca_ || chisq.probability() < minprob_){ // hit is too far from the wire or has too small a probability: disable it
        state = WireHitState::inactive;
      } else if(absdoca > mindoca_ && absdoca < maxddoca_){  // in the sweet spot: use the DOCA to sign the ambiguity
        state = doca > 0.0 ? WireHitState::right : WireHitState::left;
      } else { // hit too close to the wire to resolve ambiguity, or with a suspiciously large drift: just use the raw wire position and time to constrain the track
        state = WireHitState::null;
      }
    } else {
      state = WireHitState::inactive; // disable the hit
    }
    // update the hit
    kksh.setState(state);
  }
}
#endif
