//
// Bkg-based updater to disable background hits
//
#ifndef Mu2eKinKal_BkgStrawHitUpdater_hh
#define Mu2eKinKal_BkgStrawHitUpdater_hh
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include "Offline/Mu2eKinKal/inc/DriftInfo.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include <tuple>
#include <string>
#include <iostream>

namespace mu2e {
  class ComboHit;
 // Update based just on Bkg to the wire
  class BkgStrawHitUpdater {
    public:
      using BkgSHUConfig = std::tuple<std::string,float,std::string>;
      static std::string const& configDescription(); // description of the variables
      BkgStrawHitUpdater() : mva_(0), mvacut_(0.0) {}
      BkgStrawHitUpdater(BkgSHUConfig const& bkgshuconfig);
      WireHitState wireHitState(WireHitState const& input, KinKal::ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const;
    private:
      MVATools* mva_; // neural net calculator
      double mvacut_; // cut value to decide if drift information is usable
      WHSMask freeze_; // states to freeze
  };
}
#endif
