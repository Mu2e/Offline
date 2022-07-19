//
// enum of StrawHitUpdaters.  This needs to be extended when new updaters are added
//
#ifndef Mu2eKinKal_StrawHitUpdaters_hh
#define Mu2eKinKal_StrawHitUpdaters_hh

namespace mu2e {
  // types of updaters: these need to be extended if new updaters are defined
  namespace StrawHitUpdaters {
    enum algorithm: int {none=-1, null=0, PTCA=1, Combinatoric=10 };
    // specify which updaters operate directly on StrawHits vs StrawHitClusters
    bool updateStrawHitClusters(algorithm algo);
    bool updateStrawHits(algorithm algo);
  }
}
#endif
