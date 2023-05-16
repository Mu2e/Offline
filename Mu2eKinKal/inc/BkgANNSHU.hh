//
// Bkg-based updater to disable background hits
//
#ifndef Mu2eKinKal_BkgANNSHU_hh
#define Mu2eKinKal_BkgANNSHU_hh
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include "Offline/TrackerConditions/inc/DriftInfo.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include <tuple>
#include <string>
#include <iostream>
#include <cstddef>
#include <memory>

namespace TMVA_SOFIE_TrainBkg {
  class Session;
}

namespace mu2e {
  class ComboHit;
  // Update based just on Bkg to the wire
  class BkgANNSHU {
    public:
      using Config = std::tuple<std::string,float,std::string,int>;
      static std::string const& configDescription(); // description of the variables
      BkgANNSHU(Config const& config);
      WireHitState wireHitState(WireHitState const& input, KinKal::ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const;
    private:
      std::shared_ptr<TMVA_SOFIE_TrainBkg::Session> mva_;
      double mvacut_ =0; // cut value to decide if drift information is usable
      WHSMask freeze_; // states to freeze
      int diag_ =0; // diag print level
  };
}
#endif
