#ifndef TrkComp_HH
#define TrkComp_HH
// root 
#include "Rtypes.h"
class KalRep;
//
// Functor Class to compare tracks and find overlaps (2 tracks sharing a substantial number of hits)
// Dave Brown, LBNL 7/8/2016
namespace mu2e 
{
// struct to define overlap between 2 tracks
  struct TrkOverlap {
    UInt_t _nover; // # of active hits shared between 2 tracks
    Float_t _f1; // fraction of shared active hits in first track
    Float_t _f2; // fraction of shared active hits in second track
  };

  class TrkComp {
  public:
  // compute the number of active TrStrawHits using the same StrawHit in common between 2 tracks.  
    unsigned nOverlap(const KalRep* k1, const KalRep* k2);
    void fillOverlap(const KalRep* k1, const KalRep* k2, TrkOverlap&);
  // find all non-trivial overlaps between a collection of tracks
//    void fillOverlaps(KalRepPtrCollection const& trks,std::vector<TrkOverlap>& overlaps);
  };
}
#endif
