#include "Offline/TrkReco/inc/SimpleKalSeedSelector.hh"

namespace mu2e {

  bool SimpleKalSeedSelector::select(KalSeed const& kseed) const {
    // evaluate the momentum at t0
    if(kseed.intersections().size() > 0){
      auto const& kinter = kseed.intersections().front();
      auto mom = kinter.mom();
      auto fcon = kseed.fitConsistency();
      return mom >= minmom_ && mom <= maxmom_ && fcon >= minfcon_;
    } else
      return false;
  }

  bool SimpleKalSeedSelector::isBetter(KalSeed const& current,KalSeed const& test) const {
    double nhitfrac = 2*(test.nHits() - current.nHits())/(test.nHits() + current.nHits());
    if(nhitfrac > minsignhit_){
      // difference is significant
      return true;
    } else {
      // fall back to fit consistency
      return test.fitConsistency() > current.fitConsistency();
    }
  }


}
