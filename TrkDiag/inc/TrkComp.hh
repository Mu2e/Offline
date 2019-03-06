#ifndef TrkComp_HH
#define TrkComp_HH
// root 
#include "Rtypes.h"
class KalRep;
#include "RecoDataProducts/inc/KalSeed.hh"
//
// Functor Class to compare tracks and find overlaps (2 tracks sharing a substantial number of hits)
// Dave Brown, LBNL 7/8/2016
namespace mu2e {
  class TrkComp {
  public:
  // compute the number of active TrStrawHits using the same StrawHit in common between 2 tracks.  
    unsigned nOverlap(const KalRep* k1, const KalRep* k2);
    unsigned nOverlap(const KalSeed& k1, const KalSeed& k2);
//    void fillOverlaps(KalRepPtrCollection const& trks,std::vector<TrkOverlap>& overlaps);
  };
}
#endif
