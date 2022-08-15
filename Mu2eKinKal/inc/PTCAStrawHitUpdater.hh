//
// Simple updater of StrawHits based on Position and Time Of Closest Approach (PTCA)
// Currently this defines useable hits using a box cut, an MVA would be an improvement TODO
//
#ifndef Mu2eKinKal_PTCAStrawHitUpdater_hh
#define Mu2eKinKal_PTCAStrawHitUpdater_hh
#include "Offline/Mu2eKinKal/inc/StrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/NullHitInfo.hh"
#include <tuple>
#include <iostream>

namespace mu2e {
  // Update based just on PTCA to the wire
  class PTCAStrawHitUpdater : public StrawHitUpdater {
    public:
      using PSHUConfig = std::tuple<float,float,float,float,bool,int>;
      PTCAStrawHitUpdater() : mindoca_(0), maxdoca_(0), mindt_(0), maxdt_(0), uptca_(false), nhtmode_(NullHitInfo::none) {}
      PTCAStrawHitUpdater(PSHUConfig const& ptcashuconfig) {
        mindoca_ = std::get<0>(ptcashuconfig);
        maxdoca_ = std::get<1>(ptcashuconfig);
        mindt_ = std::get<2>(ptcashuconfig);
        maxdt_ = std::get<3>(ptcashuconfig);
        uptca_ = std::get<4>(ptcashuconfig);
        nhtmode_ = static_cast<NullHitInfo::nullTimeMode>(std::get<5>(ptcashuconfig));
        std::cout << "PTCAStrawHitUpdater " << mindoca_ << " " << maxdoca_ << " " << mindt_ << " " << maxdt_ << " " << uptca_ << " " << nhtmode_ << std::endl;
      }
      // set the state based on the current PTCA value
      WireHitState wireHitState(ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const override;
      // unassigned hit properties
      StrawHitUpdaters::algorithm algorithm() const override{ return StrawHitUpdaters::PTCA; }
      // accessors
      auto minDOCA() const { return mindoca_; }
      auto maxDOCA() const { return maxdoca_; }
      auto minDt() const { return maxdt_; }
      auto maxDt() const { return maxdt_; }
      // use biased or unbiased PTCA estimate
      bool useUnbiasedClosestApproach() const override { return uptca_; }
    private:
      double mindoca_; // minimum DOCA to use drift information
      double maxdoca_; // maximum DOCA to still use a hit; beyond this it is forced inactive
      double mindt_; // maximum dt to use drift information
      double maxdt_; // maximum dt to use drift information
      bool uptca_; // use unbiased DOCA info
      NullHitInfo::nullTimeMode nhtmode_; // use time constraint in null hits
  };
}
#endif
