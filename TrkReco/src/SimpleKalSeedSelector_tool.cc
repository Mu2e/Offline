#include "Offline/TrkReco/inc/SimpleKalSeedSelector.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"

namespace mu2e {

  bool SimpleKalSeedSelector::select(KalSeed const& kseed) const {
    // evaluate the momentum at t0
    if(kseed.intersections().size() > 0){
      auto const& kinter = kseed.intersections().front();
      auto mom = kinter.mom();
      auto fcon = kseed.fitConsistency();
      unsigned nactive = kseed.nHits(true);
      return mom >= minmom_ && mom <= maxmom_ && fcon >= minfcon_ && nactive >= minnactive_;
    } else
      return false;
  }

  bool SimpleKalSeedSelector::isBetter(KalSeed const& current,KalSeed const& test) const {
    unsigned ncurrent = current.nHits(true);
    unsigned ntest = test.nHits(true);
    if(ntest + ncurrent == 0) return test.fitConsistency() > current.fitConsistency();
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
