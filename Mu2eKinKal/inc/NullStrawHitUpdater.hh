//
// Simple class to update StrawHits as to have null ambiguity, and remove outliers
//
#ifndef Mu2eKinKal_NullStrawHitUpdater_hh
#define Mu2eKinKal_NullStrawHitUpdater_hh

#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"

namespace mu2e {
  using KinKal::ClosestApproachData;
  class ComboHit;
  class StrawResponse;
  // always set the wire hit state to null; used for seed fitting
  class NullStrawHitUpdater {
    public:
      NullStrawHitUpdater(double maxdoca) : maxdoca_(maxdoca) {}
      WireHitState wireHitState(ClosestApproachData const& tpdata ) const;
      void timeResid(ClosestApproachData const& tpdata, ComboHit const& chit,double& dt, double& dtvar) const;
      auto maxDOCA() const { return maxdoca_; }
      auto distVariance() const { return dvar_; }
      auto timeOffset() const { return toff_; }
      auto timeVariance() const { return tvar_; }
      StrawHitUpdaters::algorithm algorithm() const { return StrawHitUpdaters::null; }
   private:
      double maxdoca_; // maximum DOCA to still use a hit
      double dvar_; // variance to assign to distance
      double tvar_; // variance on time residual; should come from ComboHit FIXME
      double toff_; // time offset; this shouldn't be necessary if ComboHit were perfectly calibrated FIXME
  };
}
#endif
