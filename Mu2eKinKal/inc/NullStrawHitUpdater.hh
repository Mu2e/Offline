//
// Simple class to update StrawHits as to have null ambiguity, and remove outliers
//
#ifndef Mu2eKinKal_NullStrawHitUpdater_hh
#define Mu2eKinKal_NullStrawHitUpdater_hh
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include <tuple>

namespace mu2e {
  using KinKal::ClosestApproachData;
  class ComboHit;
  class StrawResponse;
  // always set the wire hit state to null; used for seed fitting
  class NullStrawHitUpdater {
    public:
      using NSHUConfig = std::tuple<float>;
      NullStrawHitUpdater(NSHUConfig const& nsuconfig) {
        maxdoca_ = std::get<0>(nsuconfig);
      }
      WireHitState wireHitState(WireHitState const& input,ClosestApproachData const& tpdata) const;
    private:
      double maxdoca_; // maximum DOCA to still use a hit
  };
}
#endif
