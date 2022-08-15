//
// Base class for StrawHit updaters
//
#ifndef Mu2eKinKal_StrawHitUpdater_hh
#define Mu2eKinKal_StrawHitUpdater_hh

#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/DriftInfo.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include <tuple>

namespace mu2e {
  using KinKal::ClosestApproachData;
  class StrawResponse;
  class Straw;
  class ComboHit;
  class StrawHitUpdater {
    public:
      virtual WireHitState wireHitState(ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const =0;
      virtual StrawHitUpdaters::algorithm algorithm() const =0;
      virtual bool useUnbiasedClosestApproach() const =0;
  };
}
#endif
