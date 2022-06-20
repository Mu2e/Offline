//
// Simple updater of StrawHits based on (unbiased) Distance Of Closest Approach (DOCA)
//
#ifndef Mu2eKinKal_DOCAStrawHitUpdater_hh
#define Mu2eKinKal_DOCAStrawHitUpdater_hh
#include "Offline/Mu2eKinKal/inc/WireHitStructs.hh"

namespace mu2e {
// Update based just on (unbiased) DOCA to the wire, not including this hit
  class DOCAStrawHitUpdater {
    public:
      DOCAStrawHitUpdater() : maxdoca_(-1.0), mindoca_(0), maxddoca_(-1) {}
      DOCAStrawHitUpdater(double maxdoca, double minddoca, double maxddoca) :
        maxdoca_(maxdoca), mindoca_(minddoca), maxddoca_(maxddoca) {}
     // set the state based on the current DOCA value
      WireHitState wireHitState(double doca) const;
      auto maxDOCA() const { return maxdoca_; }
      auto minDOCA() const { return mindoca_; }
      auto maxDriftDOCA() const { return maxddoca_; }
    private:
      double maxdoca_; // maximum DOCA to still use a hit
      double mindoca_; // minimum DOCA to use drift information
      double maxddoca_; // maximum DOCA to use drift information
  };

}
#endif
