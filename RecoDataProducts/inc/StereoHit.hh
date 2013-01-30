#ifndef RecoDataProducts_StereoHit_hh
#define RecoDataProducts_StereoHit_hh
//
// simple reconstruction of 2 straw hits in different views giving 3-d information from stereo 
//
// $Author: brownd $
// $Date: 2013/01/30 20:03:35 $
//
// Original author David Brown
//

// C++ includes
#include <iostream>

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
// Mu2e includes
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "TTrackerGeom/inc/TTracker.hh"

namespace mu2e {
  class StereoHit {
    public:
    // required for persistence
      StereoHit();
// construct from specified straw hits.  The collection and tracker are required
      StereoHit(StrawHitCollection const& strawhits,Tracker const& tracker, size_t ind1, size_t ind2);   // accessors
      StrawHit const& sh1(StrawHitCollection const& strawhits) const { return strawhits[_hind1]; }
      StrawHit const& sh2(StrawHitCollection const& strawhits) const { return strawhits[_hind2]; }
      Straw const& s1(StrawHitCollection const& strawhits,Tracker const& tracker) const 
	{ return tracker.getStraw(sh1(strawhits).strawIndex()); }
      Straw const& s2(StrawHitCollection const& strawhits,Tracker const& tracker) const 
	{ return tracker.getStraw(sh2(strawhits).strawIndex()); }
      CLHEP::Hep3Vector const& pos() const { return _pos; }
      float dist() const { return _dist; }
      float time() const { return _time; }
      float dt() const { return _dt; } // signed t2 -t1
      float wdist1() const { return _wd1; }
      float wdist2() const { return _wd2; }
      size_t hitIndex1() const { return _hind1; }
      size_t hitIndex2() const { return _hind2; }
    private:
      size_t _hind1, _hind2; // indices into the straw hit container for the 2 hits making up this stereo hit
      CLHEP::Hep3Vector _pos; // position in tracker coordinates
      float _dist; // transverse separation between the 2 hit wires at their POCA
      float _time; // average time of the 2 hits
      float _dt; // time separation of the 2 hits
      float _wd1, _wd2; // distance along the wire of the POCA for the two hits
    };
}
#endif

