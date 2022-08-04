//
// Simple updater of StrawHits based on Position and Time Of Closest Approach (PTCA)
// Currently this defines useable hits using a box cut, an MVA would be an improvement TODO
//
#ifndef Mu2eKinKal_PTCAStrawHitUpdater_hh
#define Mu2eKinKal_PTCAStrawHitUpdater_hh
#include "Offline/Mu2eKinKal/inc/StrawHitUpdater.hh"
#include <tuple>
#include <iostream>

namespace mu2e {
  // Update based just on PTCA to the wire
  class PTCAStrawHitUpdater : public StrawHitUpdater {
    public:
      enum nullTimeMode{none=0,usedoca,usecombo};
      using DSHUConfig = std::tuple<float,float,float,float,bool,int>;
      PTCAStrawHitUpdater() : mindoca_(0), maxdoca_(0), mindt_(0), maxdt_(0), uptca_(false), nulltime_(none), dvar_(0) {}
      PTCAStrawHitUpdater(DSHUConfig const& dshusetting) {
        mindoca_ = std::get<0>(dshusetting);
        maxdoca_ = std::get<1>(dshusetting);
        mindt_ = std::get<2>(dshusetting);
        maxdt_ = std::get<3>(dshusetting);
        uptca_ = std::get<4>(dshusetting);
        nulltime_ = static_cast<nullTimeMode>(std::get<5>(dshusetting));
        static double invthree(1.0/3.0);
        dvar_ = invthree*mindoca_*mindoca_;
        std::cout << "PTCAStrawHitUpdater " << mindoca_ << " " << maxdoca_ << " " << mindt_ << " " << maxdt_ << " " << uptca_ << " " << nulltime_ << std::endl;
      }
      // set the state based on the current PTCA value
      WireHitState wireHitState(ClosestApproachData const& tpdata, Straw const& straw, StrawResponse const& sresponse) const override;
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
      nullTimeMode nulltime_; // use time constraint in null hits
      double dvar_; // null hit distance variance
  };
}
#endif
