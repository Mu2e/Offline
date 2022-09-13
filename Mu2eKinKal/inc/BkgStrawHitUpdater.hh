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
      using BkgSHUConfig = std::tuple<std::string,float,std::string,int>;
      static std::string const& configDescription(); // description of the variables
      BkgStrawHitUpdater(BkgSHUConfig const& bkgshuconfig);
      BkgStrawHitUpdater(BkgStrawHitUpdater const& other) :  mvacut_(other.mvacut_), freeze_(other.freeze_), diag_(other.diag_) {
        if(other.mva_) mva_ = new MVATools(*other.mva_);
      }
      ~BkgStrawHitUpdater() { delete mva_; }
      WireHitState wireHitState(WireHitState const& input, KinKal::ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const;
    private:
      MVATools* mva_ =0; // neural net calculator
      double mvacut_ =0; // cut value to decide if drift information is usable
      WHSMask freeze_; // states to freeze
      int diag_ =0; // diag print level
  };
}
#endif
