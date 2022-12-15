//
// Bkg-based updater to disable background hits
//
#ifndef Mu2eKinKal_BkgANNSHU_hh
#define Mu2eKinKal_BkgANNSHU_hh
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
 // Update based just on Bkg to the wire
  class BkgANNSHU {
    public:
      using Config = std::tuple<std::string,float,std::string,int>;
      static std::string const& configDescription(); // description of the variables
      BkgANNSHU(Config const& config);
      BkgANNSHU(BkgANNSHU const& other) :  mvacut_(other.mvacut_), freeze_(other.freeze_), diag_(other.diag_) {
        if(other.mva_) mva_ = new MVATools(*other.mva_);
      }
      ~BkgANNSHU() { delete mva_; }
      WireHitState wireHitState(WireHitState const& input, KinKal::ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const;
    private:
      MVATools* mva_ =nullptr; // neural net calculator
      double mvacut_ =0; // cut value to decide if drift information is usable
      WHSMask freeze_; // states to freeze
      int diag_ =0; // diag print level
  };
}
#endif
