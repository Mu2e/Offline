//
// Simple class to update StrawHits as to have null ambiguity, and remove outliers
//
#ifndef Mu2eKinKal_NullStrawHitUpdater_hh
#define Mu2eKinKal_NullStrawHitUpdater_hh

#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"

namespace mu2e {
  // always set the wire hit state to null; used for seed fitting
  class NullStrawHitUpdater {
    public:
      NullStrawHitUpdater(double maxdoca) : maxdoca_(maxdoca) {}
      WireHitState wireHitState(KinKal::ClosestApproachData const& tpdata ) const;
      auto maxDOCA() const { return maxdoca_; }
    private:
      double maxdoca_; // maximum DOCA to still use a hit
  };
}
#endif
