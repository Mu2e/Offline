#ifndef RecoDataProducts_LineSeed_hh
#define RecoDataProducts_LineSeed_hh
//
// Seed for a straight line track
// Richie Bonventre
//

// Mu2e includes
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <vector>

namespace mu2e {
  struct LineSeed {

    CLHEP::Hep3Vector _seedDir;
    CLHEP::Hep3Vector _seedInt; // y intercept
    double _t0;
    int _seedSize;
    bool _converged;
    art::Ptr<TimeCluster>    _timeCluster; // associated time cluster
    std::vector<StrawHitIndex> _strawHitIdxs; // associated straw hits: can be empty
  };
   typedef std::vector<mu2e::LineSeed> LineSeedCollection;
} // namespace mu2e

#endif /* RecoDataProducts_LneSeed_hh */
