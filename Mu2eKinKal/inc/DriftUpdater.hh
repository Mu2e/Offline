//
// Interface for updaters that work on StrawHit drift
//
#ifndef Mu2eKinKal_DriftUpdater_hh
#define Mu2eKinKal_Driftpdater_hh
#include "Offline/TrackerConditions/inc/StrawResponse.hh"

namespace mu2e {
// Update based just on (unbiased) DOCA to the wire, not including this hit
  class DOCAStrawHitUpdater {
   public:
     DriftUpdater(){}

      StrawHitUpdaters::algorithm algorithm() const { return StrawHitUpdaters::DOCA; }
     // set the state based on the current DOCA value
      WireHitState wireHitState(StrawResponse const& sresponse, KinKal::ClosestApproachData const& tpdata) const;
      auto maxDOCA() const { return maxdoca_; }
      virtual double minDOCA() const=0;
      auto maxDriftDOCA() const { return maxddoca_; }
      auto maxDt() const { return maxdt_; }
  };

}
#endif
