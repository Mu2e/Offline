#ifndef Mu2eKinKal_KKCombinatoricUpdater_hh
#define Mu2eKinKal_KKCombinatoricUpdater_hh
//
//  StrawHitGroup updating using an exhaustive combinatoric algorithm, following the BTrk PanelAmbigResolver algorithm
//
#include "KinKal/Detector/WireHitStructs.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"

namespace mu2e {
  // interface for updating individual straw hits
  template <class KTRAJ> class KKStrawHitGroup;
  class KKCombinatoricUpdater {
    public:
      using PHState = std::vector<WireHitState>; // complete state of every StrawHit in a Combinatoric
     KKCombinatoricUpdater() {;}
//      KKCombinatoricUpdater(){;}
      template <class KTRAJ> void update(KKStrawHitGroup<KTRAJ>& swh) const;
    private:
  };

  template <class KTRAJ> void KKCombinatoricUpdater::update(KKStrawHitGroup<KTRAJ>& kkph) const {
    using KinKal::ClosestApproachData;
    using KinKal::WireHitState;
  }
}
#endif
