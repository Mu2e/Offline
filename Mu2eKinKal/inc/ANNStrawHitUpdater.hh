//
// ANN-based StrawHits based on Position and Time Of Closest Approach (ANN)
// Currently this defines useable hits using a box cut, an MVA would be an improvement TODO
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
      using ANNSHUConfig = std::tuple<std::string,float,int,float>;
      ANNStrawHitUpdater() : mva_(0), mvacut_(0.0), nhmode_(WireHitState::none), dvar_(0) {}
      ANNStrawHitUpdater(ANNSHUConfig const& annshuconfig) {
        mva_  = new MVATools(std::get<0>(annshuconfig));
        mvacut_ = std::get<1>(annshuconfig);
        nhmode_ = static_cast<WireHitState::NHMode>(std::get<2>(annshuconfig));
        dvar_ = std::get<3>(annshuconfig);
        std::cout << "ANNStrawHitUpdater " << " anncut " << mvacut_ << " null dvar, mode" << dvar_ << nhmode_ << std::endl;
        mva_->initMVA();
        mva_->showMVA();
      }
      // set the state based on the current ANN value
      WireHitState wireHitState(KinKal::ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const;
    private:
      MVATools* mva_; // neural net calculator
      double mvacut_; // cut value to decide if drift information is usable
      WireHitState::NHMode nhmode_; // Null hit mode
      double dvar_; // null hit distance variance
  };
}
#endif
