//
// ANN-based  ambiguity updater.  This assigns LR ambiguity based on the predicted accuracy of the DOCA sign
//
#ifndef Mu2eKinKal_ANNStrawHitUpdater_hh
#define Mu2eKinKal_ANNStrawHitUpdater_hh
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include "Offline/Mu2eKinKal/inc/DriftInfo.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include <tuple>
#include <string>
#include <iostream>
#include <cstddef>

namespace mu2e {
  class ComboHit;
  // Update based just on ANN to the wire
  class ANNStrawHitUpdater {
    public:
      using ANNSHUConfig = std::tuple<std::string,float,float,std::string,int>;
      ANNStrawHitUpdater(ANNSHUConfig const& annshuconfig);
      ANNStrawHitUpdater(ANNStrawHitUpdater const& other) : mvacut_(other.mvacut_), nulldoca_(other.nulldoca_), freeze_(other.freeze_) {
        if(other.mva_) mva_ = new MVATools(*other.mva_);
      }
      ~ANNStrawHitUpdater() { delete mva_; }
      WireHitState wireHitState(WireHitState const& input, KinKal::ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const;
      static std::string const& configDescription(); // description of the variables
    private:
      MVATools* mva_ = nullptr; // neural net calculator
      double mvacut_ =0; // cut value to decide if drift information is usable
      double nulldoca_ =2.5; // null hit doca
      WHSMask freeze_; // states to freeze
      int diag_; // diag print level
  };
}
#endif
