//
// Simple class to update StrawHits as to have null ambiguity, and remove outliers
//
#ifndef Mu2eKinKal_NullStrawHitUpdater_hh
#define Mu2eKinKal_NullStrawHitUpdater_hh

#include "KinKal/Detector/WireHitStructs.hh"

namespace mu2e {
  // always set the wire hit state to null; used for seed fitting
  class NullStrawHitUpdater {
    public:
      NullStrawHitUpdater() : maxdoca_(1.0e6) {} // default is to turn all hits null
      NullStrawHitUpdater(double maxdoca) : maxdoca_(maxdoca) {}
      KinKal::WireHitState wireHitState(double doca ) const {
        if(fabs(doca) > maxdoca_)
          return KinKal::WireHitState(KinKal::WireHitState::inactive);
        else
          return KinKal::WireHitState(KinKal::WireHitState::null);
      }
      auto maxDOCA() const { return maxdoca_; }
    private:
      double maxdoca_; // maximum DOCA to still use a hit
  };
}
#endif
