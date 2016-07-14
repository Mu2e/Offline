//
// Functor Class to compare tracks and find overlaps (2 tracks sharing a substantial number of hits)
// Dave Brown, LBNL 7/8/2016
#include "TrkDiag/inc/TrkComp.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include <functional>

namespace mu2e {

// count use of the same StrawHitis in different tracks
  unsigned TrkComp::nOverlap(const KalRep* k1, const KalRep* k2){
    unsigned nover(0);
  // count shared StrawHits.  Note this currently requires a downcast, it should
  // be supported by the TrkHit base class, FIXME!
    TrkStrawHitVector tshv1;
    convert(k1->hitVector(),tshv1);
    TrkStrawHitVector tshv2;
    convert(k2->hitVector(),tshv2);
    // sort these by StrawHit index, which is unique
    indexcomp icomp;
    std::sort(tshv1.begin(),tshv1.end(),icomp);
    std::sort(tshv2.begin(),tshv2.end(),icomp);
    // find the 1st active hit in each vector
    auto itsh1 = tshv1.begin();
    while(itsh1 != tshv1.end() && !(*itsh1)->isActive() )++itsh1;
    auto itsh2 = tshv2.begin();
    while(itsh2 != tshv2.end() && !(*itsh2)->isActive() )++itsh2;
    // loop over all potential overlapps
    while(itsh1 != tshv1.end() && itsh2 != tshv2.end()){
      // check for overlap
      if((*itsh1)->index() == (*itsh2)->index()) {
	// same straw hit: increment overlap and both iterators
	++nover; ++itsh1; ++itsh2;
      	while(itsh1 != tshv1.end() && !(*itsh1)->isActive() )++itsh1;
	while(itsh2 != tshv2.end() && !(*itsh2)->isActive() )++itsh2;
      } else if((*itsh1)->index() < (*itsh2)->index()) {
      // move to next potential overlap
	while(itsh1 != tshv1.end() && (*itsh1)->index() < (*itsh2)->index() ){
	  ++itsh1;
	  while(itsh1 != tshv1.end() && !(*itsh1)->isActive() )++itsh1;
	}
      } else {
	while(itsh2 != tshv2.end() && (*itsh2)->index() < (*itsh1)->index() ){
	  ++itsh2;
	  while(itsh2 != tshv2.end() && !(*itsh2)->isActive() )++itsh2;
	}
      }
    }
    return nover;
  }

  void TrkComp::fillOverlap(const KalRep* k1, const KalRep* k2,TrkOverlap& overlap){
  // count overlaps
    unsigned nover = nOverlap(k1,k2);
    overlap._nover = nover;
    overlap._f1 = overlap._f2 = 0;
    // create output object
    if(nover > 0){
    // count active hits in each track
      overlap._f1 = nover/(float)k1->nActive();
      overlap._f2 = nover/(float)k2->nActive();
    }
  }
}

