//
// Simple class to update StrawHits as to have null ambiguity, and remove outliers
//
#ifndef Mu2eKinKal_NullStrawHitUpdater_hh
#define Mu2eKinKal_NullStrawHitUpdater_hh

#include "Offline/Mu2eKinKal/inc/StrawHitUpdater.hh"
#include <tuple>

namespace mu2e {
  using KinKal::ClosestApproachData;
  class ComboHit;
  class StrawResponse;
  // always set the wire hit state to null; used for seed fitting
  class NullStrawHitUpdater : public StrawHitUpdater {
    public:
      using NSHUConfig = std::tuple<float,float,float,bool>;
      NullStrawHitUpdater(NSHUConfig const& nsuconfig) {
        maxdoca_ = std::get<0>(nsuconfig);
        dvar_ = std::get<1>(nsuconfig);
        tvar_ = std::get<2>(nsuconfig);
        uptca_ = std::get<3>(nsuconfig);
      }
      WireHitState wireHitState(ClosestApproachData const& tpdata, Straw const& straw, StrawResponse const& sresponse ) const override;
      StrawHitUpdaters::algorithm algorithm() const override { return StrawHitUpdaters::null; }
      bool useUnbiasedClosestApproach() const override { return uptca_; }
    private:
      double maxdoca_; // maximum DOCA to still use a hit
      double dvar_; // variance to assign to distance
      double tvar_; // time variance; should come from ComboHit
      bool uptca_; // use unbiased PTCA info
  };
}
#endif
