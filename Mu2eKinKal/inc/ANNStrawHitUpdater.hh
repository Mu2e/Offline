//
// ANN-based StrawHits based on Position and Time Of Closest Approach (ANN)
// Currently this defines useable hits using a box cut, an MVA would be an improvement TODO
//
#ifndef Mu2eKinKal_ANNStrawHitUpdater_hh
#define Mu2eKinKal_ANNStrawHitUpdater_hh
#include "Offline/Mu2eKinKal/inc/StrawHitUpdater.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/Mu2eKinKal/inc/NullHitInfo.hh"
#include <tuple>
#include <string>
#include <iostream>

namespace mu2e {
  // Update based just on ANN to the wire
  class ANNStrawHitUpdater : public StrawHitUpdater {
    public:
      using ANNSHUConfig = std::tuple<std::string,float,bool,int,float>;
      ANNStrawHitUpdater() : mva_(0), mvacut_(0.0), uptca_(false), nhtmode_(NullHitInfo::none), dvar_(0) {}
      ANNStrawHitUpdater(ANNSHUConfig const& annshuconfig) {
        mva_  = new MVATools(std::get<0>(annshuconfig));
        mvacut_ = std::get<1>(annshuconfig);
        uptca_ = std::get<2>(annshuconfig);
        nhtmode_ = static_cast<NullHitInfo::nullTimeMode>(std::get<3>(annshuconfig));
        dvar_ = std::get<4>(annshuconfig);
        std::cout << "ANNStrawHitUpdater " << " anncut " << mvacut_ << " null dvar, mode" << dvar_ << nhtmode_ << std::endl;
        mva_->initMVA();
        mva_->showMVA();
      }
      // set the state based on the current ANN value
      WireHitState wireHitState(ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const override;
      // unassigned hit properties
      StrawHitUpdaters::algorithm algorithm() const override{ return StrawHitUpdaters::ANN; }
      // accessors
      // use biased or unbiased ANN estimate
      bool useUnbiasedClosestApproach() const override { return uptca_; }
    private:
      MVATools* mva_; // neural net calculator
      double mvacut_; // cut value to decide if drift information is usable
      bool uptca_; // use unbiased DOCA info
      NullHitInfo::nullTimeMode nhtmode_; // Null hit mode
      double dvar_; // null hit distance variance
  };
}
#endif
