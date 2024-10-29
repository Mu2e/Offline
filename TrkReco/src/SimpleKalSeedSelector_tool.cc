#include "Offline/TrkReco/inc/SimpleKalSeedSelector.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"

namespace mu2e {

  bool SimpleKalSeedSelector::select(KalSeed const& kseed) const {
    // evaluate the momentum at t0
    if(kseed.intersections().size() > 0){
      auto const& kinter = kseed.intersections().front();
      auto mom = kinter.mom();
      auto fcon = kseed.fitConsistency();
      unsigned nactive =0;
      for (auto const& hit : kseed.hits()){
        if (hit.strawHitState() > WireHitState::inactive) ++nactive;
      }
      return mom >= minmom_ && mom <= maxmom_ && fcon >= minfcon_ && nactive >= minnactive_;
    } else
      return false;
  }

  bool SimpleKalSeedSelector::isBetter(KalSeed const& current,KalSeed const& test) const {
    unsigned ncurrent =0;
    for (auto const& hit : current.hits()){
      if (hit.strawHitState() > WireHitState::inactive) ++ncurrent;
    }

    unsigned ntest =0;
    for (auto const& hit : test.hits()){
      if (hit.strawHitState() > WireHitState::inactive) ++ntest;
    }

    float nhitfrac = 2*float(ntest - ncurrent)/float(ntest + ncurrent);
    if(fabs(nhitfrac) > minsignhit_){
      // hit difference is significant;
      return ntest > ncurrent;
    } else {
      // fall back to fit consistency
      return test.fitConsistency() > current.fitConsistency();
    }
  }


}
