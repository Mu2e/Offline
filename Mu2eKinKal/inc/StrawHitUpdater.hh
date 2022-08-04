//
// Base class for StrawHit updaters
//
#ifndef Mu2eKinKal_StrawHitUpdater_hh
#define Mu2eKinKal_StrawHitUpdater_hh

#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include <tuple>

namespace mu2e {
  using KinKal::ClosestApproachData;
  class StrawResponse;
  class Straw;
  // always set the wire hit state to null; used for seed fitting
  class StrawHitUpdater {
    public:
      virtual WireHitState wireHitState(ClosestApproachData const& tpdata, Straw const& straw, StrawResponse const& sresponse ) const =0;;
      virtual StrawHitUpdaters::algorithm algorithm() const =0;
      virtual bool useUnbiasedClosestApproach() const =0;
  };
}
#endif
