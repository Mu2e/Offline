//
// Simple updater of StrawHits based on box cuts of Closest Approach (CA) information.  This mimics BTrk
//
#ifndef Mu2eKinKal_CAStrawHitUpdater_hh
#define Mu2eKinKal_CAStrawHitUpdater_hh
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include <tuple>
#include <string>
#include <iostream>

namespace mu2e {
  // Update based just on PTCA to the wire
  class CAStrawHitUpdater {
    public:
      using CASHUConfig = std::tuple<float,float,float,float,std::string>;
      CAStrawHitUpdater() : mindoca_(0), maxdoca_(0), mindt_(0), maxdt_(0) {}
      CAStrawHitUpdater(CASHUConfig const& cashuconfig);
      // set the state based on the current PTCA value
      WireHitState wireHitState(WireHitState const& input, KinKal::ClosestApproachData const& tpdata) const;
      // accessors
      auto minDOCA() const { return mindoca_; }
      auto maxDOCA() const { return maxdoca_; }
      auto minDt() const { return maxdt_; }
      auto maxDt() const { return maxdt_; }
    private:
      double mindoca_; // minimum DOCA to use drift information
      double maxdoca_; // maximum DOCA to still use a hit; beyond this it is forced inactive
      double mindt_; // maximum dt to use drift information
      double maxdt_; // maximum dt to use drift information
      WHSMask freeze_; // states to freeze
  };
}
#endif
