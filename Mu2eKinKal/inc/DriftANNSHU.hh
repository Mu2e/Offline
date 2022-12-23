//
// ANN-based  ambiguity updater.  This assigns LR ambiguity based on the predicted accuracy of the DOCA sign
//
#ifndef Mu2eKinKal_DriftANNSHU_hh
#define Mu2eKinKal_DriftANNSHU_hh
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include "Offline/TrackerConditions/inc/DriftInfo.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include <tuple>
#include <string>
#include <iostream>
#include <memory>
#include <cstddef>

namespace mu2e {
  class ComboHit;
  // Update based just on ANN to the wire
  class DriftANNSHU {
    public:
      using Config = std::tuple<std::string,float,std::string,std::string,std::string,std::string,int>;
      DriftANNSHU(Config const& config);
      WireHitState wireHitState(WireHitState const& input, KinKal::ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const;
      static std::string const& configDescription(); // description of the variables
    private:
      std::shared_ptr<MVATools> mva_;
      double mvacut_ =0; // cut value to decide if drift information is usable
      WireHitState::NullDistVar nulldvar_; // null hit doca
      WireHitState::TOTUse totuse_ = WireHitState::all; // use TOT time as a residual for all hits
      WHSMask allowed_; // allowed states
      WHSMask freeze_; // states to freeze
      int diag_; // diag print level
  };
}
#endif
