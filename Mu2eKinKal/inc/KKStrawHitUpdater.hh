#ifndef Mu2eKinKal_KKStrawHitUpdater_hh
#define Mu2eKinKal_KKStrawHitUpdater_hh

#include "KinKal/Detector/WireHitStructs.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include <limits>
using KinKal::ClosestApproachData;
  using KinKal::WireHitState;

namespace mu2e {
// interface for updating straw hits
  class KKStrawHitUpdater {
    public:
    virtual void updateState(WireHitState& hitstate, ClosestApproachData const& poca) const = 0;
  };
  // simple implementation of the above
  class KKSimpleStrawHitUpdater : public KKStrawHitUpdater {
    public:
      KKSimpleStrawHitUpdater() : mindoca_(std::numeric_limits<float>::max()), maxdoca_(-1.0), nulldim_(WireHitState::both) {}
      KKSimpleStrawHitUpdater(double mindoca, double maxdoca, WireHitState::Dimension nulldim) : mindoca_(mindoca), maxdoca_(maxdoca), nulldim_(nulldim) {}
      void updateState(WireHitState& hitstate, ClosestApproachData const& poca) const override;
    private:
      double mindoca_; // minimum DOCA value to use drift information
      double maxdoca_; // maximum DOCA to still use a hit
      WireHitState::Dimension nulldim_; // constrain dimension for null hits
  };
}
#endif
