#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
namespace mu2e {
  namespace StrawHitUpdaters {
    bool updateStrawHitClusters(algorithm algo) { return algo >= Combinatoric; }
    bool updateStrawHits(algorithm algo) { return algo > none && algo < Combinatoric; }
  }
}
