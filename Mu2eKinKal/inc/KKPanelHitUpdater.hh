#ifndef Mu2eKinKal_KKPanelHitUpdater_hh
#define Mu2eKinKal_KKPanelHitUpdater_hh
//
//  Panel hit updating, following the BTrk PanelAmbigResolver algorithm
//
#include "KinKal/Detector/WireHitStructs.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"

namespace mu2e {
  // interface for updating individual straw hits
  template <class KTRAJ> class KKPanelHit;
  class KKPanelHitUpdater {
    public:
      using PHState = std::array<WireHitState,MAXNHIT>; // complete state of every StrawHit in a PanelHit
     KKPanelHitUpdater() {;}
//      KKPanelHitUpdater(){;}
      template <class KTRAJ> void update(KKPanelHit<KTRAJ>& swh) const;
    private:
  };

  template <class KTRAJ> void KKPanelHitUpdater::update(KKPanelHit<KTRAJ>& kkph) const {
    using KinKal::ClosestApproachData;
    using KinKal::WireHitState;
  }
}
#endif
