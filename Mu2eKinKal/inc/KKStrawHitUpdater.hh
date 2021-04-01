#ifndef Mu2eKinKal_KKStrawHitUpdater_hh
#define Mu2eKinKal_KKStrawHitUpdater_hh

#include "KinKal/Detector/WireHitStructs.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include <limits>
using KinKal::ClosestApproachData;
using KinKal::WireHitState;

namespace mu2e {
  // class for updating straw hits
  class KKStrawHitUpdater {
    public:
      KKStrawHitUpdater() : mindoca_(std::numeric_limits<float>::max()), maxdoca_(-1.0), nulldim_(WireHitState::both) {}
      KKStrawHitUpdater(double mindoca, double maxdoca, WireHitState::Dimension nulldim) : mindoca_(mindoca), maxdoca_(maxdoca), nulldim_(nulldim) {}
      void updateState(WireHitState& hitstate, ClosestApproachData const& poca) const;
    private:
      double mindoca_; // minimum DOCA value to use drift information
      double maxdoca_; // maximum DOCA to still use a hit
      WireHitState::Dimension nulldim_; // constrain dimension for null hits
  };
}
#endif
