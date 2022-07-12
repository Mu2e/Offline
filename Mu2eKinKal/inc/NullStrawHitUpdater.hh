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
      using NSHUConfig = std::tuple<float,float>;
      NullStrawHitUpdater(NSHUConfig const& nsuconfig) {
        maxdoca_ = std::get<0>(nsuconfig);
        dvar_ = std::get<1>(nsuconfig);
      }
      WireHitState wireHitState(ClosestApproachData const& tpdata ) const;
      void timeResid(ClosestApproachData const& tpdata, ComboHit const& chit,double& dt, double& dtvar) const;
      auto maxDOCA() const { return maxdoca_; }
      auto distVariance() const { return dvar_; }
      StrawHitUpdaters::algorithm algorithm() const { return StrawHitUpdaters::null; }
   private:
      double maxdoca_; // maximum DOCA to still use a hit
      double dvar_; // variance to assign to distance
  };
}
#endif
