//
// Bkg-based updater to disable background hits
//
#ifndef Mu2eKinKal_BkgStrawHitUpdater_hh
#define Mu2eKinKal_BkgStrawHitUpdater_hh
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
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
      using BkgSHUConfig = std::tuple<std::string,float,bool>;
      BkgStrawHitUpdater() : mva_(0), mvacut_(0.0), freeze_(false) {}
      BkgStrawHitUpdater(BkgSHUConfig const& bkgshuconfig) {
        mva_  = new MVATools(std::get<0>(bkgshuconfig));
        mvacut_ = std::get<1>(bkgshuconfig);
        freeze_ = std::get<2>(bkgshuconfig);
        std::cout << "BkgStrawHitUpdater " << " bkgcut " << mvacut_ << " freeze " << freeze_ << std::endl;
        mva_->initMVA();
        mva_->showMVA();
      }
      WireHitState wireHitState(WireHitState const& input, KinKal::ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const;
    private:
      MVATools* mva_; // neural net calculator
      double mvacut_; // cut value to decide if drift information is usable
      bool freeze_; // freeze disabled hits or not
  };
}
#endif
