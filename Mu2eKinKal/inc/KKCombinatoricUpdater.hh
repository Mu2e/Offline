#ifndef Mu2eKinKal_KKCombinatoricUpdater_hh
#define Mu2eKinKal_KKCombinatoricUpdater_hh
//
//  StrawHitSet updating using an exhaustive combinatoric algorithm, following the BTrk PanelAmbigResolver algorithm
//
#include "KinKal/Detector/WireHitStructs.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"

namespace mu2e {
  // interface for updating individual straw hits
  template <class KTRAJ> class KKStrawHitSet;
  class KKCombinatoricUpdater {
    public:
      using PHState = std::vector<WireHitState>; // complete state of every StrawHit in a Combinatoric
     KKCombinatoricUpdater() {;}
//      KKCombinatoricUpdater(){;}
      template <class KTRAJ> void update(KKStrawHitSet<KTRAJ>& swh) const;
    private:
  };

  template <class KTRAJ> void KKCombinatoricUpdater::update(KKStrawHitSet<KTRAJ>& kkph) const {
    using KinKal::ClosestApproachData;
    using KinKal::WireHitState;
  }
}
#endif
