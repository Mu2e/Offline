#ifndef Mu2eKinKal_CADSHU_hh
#define Mu2eKinKal_CADSHU_hh
//
// Simple updater of StrawHits based on box cuts of Closest Approach (CA) and drift information
//
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include "Offline/TrackerConditions/inc/DriftInfo.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include <tuple>
#include <string>
#include <iostream>

namespace mu2e {
  // Update based just on PTCA to the wire
  class CADSHU {
    public:
      using Config = std::tuple<float,float,float,float,std::string,std::string,std::string,std::string,int>;
      CADSHU(Config const& config);
      static std::string const& configDescription(); // description of the variables
      // set the state based on the current PTCA value
      WireHitState wireHitState(WireHitState const& input, KinKal::ClosestApproachData const& tpdata,DriftInfo const& dinfo) const;
    private:
      double maxdoca_ =0; // maximum DOCA to use hit
      double maxdvar_ =0; // maximum DOCA variance to use hit
      double minrdrift_ =0; // minimum rdrift to use drift information
      double maxrdrift_ =0; // maximum rdrift to use hit
      WHSMask allowed_; // allowed states
      WHSMask freeze_; // states to freeze
      WireHitState::NullDistVar nulldvar_; // null hit distance variance setting
      WireHitState::TOTUse totuse_; // TOT time constraint use
      int diag_ =0; // diag print level
  };
}
#endif
