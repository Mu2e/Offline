//
// Simple updater of StrawHits based on Position and Time Of Closest Approach (PTCA)
// Currently this defines useable hits using a box cut, an MVA would be an improvement TODO
//
#ifndef Mu2eKinKal_PTCAStrawHitUpdater_hh
#define Mu2eKinKal_PTCAStrawHitUpdater_hh
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include <tuple>
#include <iostream>

namespace mu2e {
  // Update based just on PTCA to the wire
  class PTCAStrawHitUpdater {
    public:
      using PSHUConfig = std::tuple<float,float,float,float,int>;
      PTCAStrawHitUpdater() : mindoca_(0), maxdoca_(0), mindt_(0), maxdt_(0), nhmode_(WireHitState::combo) {}
      PTCAStrawHitUpdater(PSHUConfig const& ptcashuconfig) {
        mindoca_ = std::get<0>(ptcashuconfig);
        maxdoca_ = std::get<1>(ptcashuconfig);
        mindt_ = std::get<2>(ptcashuconfig);
        maxdt_ = std::get<3>(ptcashuconfig);
        nhmode_ = static_cast<WireHitState::NHMode>(std::get<4>(ptcashuconfig));
        dvar_ = mindoca_*mindoca_/3.0;
       std::cout << "PTCAStrawHitUpdater " << mindoca_ << " " << maxdoca_ << " " << mindt_ << " " << maxdt_ << " " << nhmode_ << std::endl;
      }
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
      double dvar_; // distance variance
      WireHitState::NHMode nhmode_; // null hit mode
  };
}
#endif
