//
// Simple updater of StrawHits based on (unbiased) Distance Of Closest Approach (DOCA)
//
#ifndef Mu2eKinKal_DOCAStrawHitUpdater_hh
#define Mu2eKinKal_DOCAStrawHitUpdater_hh
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include <tuple>

namespace mu2e {
// Update based just on (unbiased) DOCA to the wire, not including this hit
  class DOCAStrawHitUpdater {
   public:
     using DSHUConfig = std::tuple<float,float,float,bool>;
     DOCAStrawHitUpdater() : maxdoca_(-1.0), mindoca_(0), maxdt_(-1), uptca_(true) {}
     DOCAStrawHitUpdater(DSHUConfig const& dshusetting) {
       maxdoca_ = std::get<0>(dshusetting);
       mindoca_ = std::get<1>(dshusetting);
       maxdt_ = std::get<2>(dshusetting);
       uptca_ = std::get<3>(dshusetting);
       static double invthree(1.0/3.0);
       dvar_ = invthree*mindoca_*mindoca_;
     }
     // set the state based on the current DOCA value
      WireHitState wireHitState(StrawResponse const& sresponse, KinKal::ClosestApproachData const& tpdata) const;
      auto maxDOCA() const { return maxdoca_; }
      auto minDOCA() const { return mindoca_; }
      auto distVariance() const { return dvar_; }
      auto maxDt() const { return maxdt_; }
      auto useUnbiasedDOCA() const { return uptca_; }
      StrawHitUpdaters::algorithm algorithm() const { return StrawHitUpdaters::DOCA; }
   private:
      double maxdoca_; // maximum DOCA to still use a hit; beyond this it is forced inactive
      double mindoca_; // minimum DOCA to use drift information
      double maxdt_; // maximum dt to use drift information
      bool uptca_; // use unbiased DOCA info
      double dvar_; // distance variance
  };
}
#endif
