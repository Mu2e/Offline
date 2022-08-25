//
// ANN-based  ambiguity updater.  This assigns LR ambiguity based on the predicted accuracy of the DOCA sign
//
#ifndef Mu2eKinKal_ANNStrawHitUpdater_hh
#define Mu2eKinKal_ANNStrawHitUpdater_hh
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
 // Update based just on ANN to the wire
  class ANNStrawHitUpdater {
    public:
      using ANNSHUConfig = std::tuple<std::string,float,int,float,bool>;
      ANNStrawHitUpdater() : mva_(0), mvacut_(0.0), nhmode_(WireHitState::none), dvar_(0), freeze_(false) {}
      ANNStrawHitUpdater(ANNSHUConfig const& annshuconfig) {
        mva_  = new MVATools(std::get<0>(annshuconfig));
        mvacut_ = std::get<1>(annshuconfig);
        nhmode_ = static_cast<WireHitState::NHMode>(std::get<2>(annshuconfig));
        dvar_ = std::get<3>(annshuconfig);
        freeze_ = std::get<4>(annshuconfig);
        std::cout << "ANNStrawHitUpdater " << " anncut " << mvacut_ << " null dvar, mode" << dvar_ << nhmode_ << " freeze " << freeze_ << std::endl;
        mva_->initMVA();
        mva_->showMVA();
      }
      WireHitState wireHitState(WireHitState const& input, KinKal::ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const;
    private:
      MVATools* mva_; // neural net calculator
      double mvacut_; // cut value to decide if drift information is usable
      WireHitState::NHMode nhmode_; // Null hit mode
      double dvar_; // null hit distance variance
      bool freeze_; // freeze drift states
  };
}
#endif
