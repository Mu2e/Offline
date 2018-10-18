#ifndef RecoDataProducts_StereoHit_hh
#define RecoDataProducts_StereoHit_hh
#include "RecoDataProducts/inc/ComboHit.hh"
namespace mu2e {
  typedef ComboHit StereoHit;
  typedef ComboHitCollection StereoHitCollection;
}
#endif
////
//// simple reconstruction of 2 straw hits in different views giving 3-d information from stereo
////
//// $Author: brownd $
//// $Date: 2013/03/08 04:29:49 $
////
//// Original author David Brown
////
//
//// C++ includes
//#include <iostream>
//
//// CLHEP includes
//#include "CLHEP/Vector/ThreeVector.h"
//// Mu2e includes
//#include "RecoDataProducts/inc/StrawHit.hh"
//#include "RecoDataProducts/inc/StrawHitCollection.hh"
//#include "DataProducts/inc/PanelId.hh"
//// C++
//#include <vector>
//
//namespace mu2e {
//  class StereoHit {
//    public:
//    // required for persistence
//      StereoHit();
//// construct from specified straw hits.  The collection and tracker are required
//      StereoHit(StrawHitCollection const& strawhits,Tracker const& tracker, std::size_t ind1, std::size_t ind2);
//// accessors
//      StrawHit const& sh1(StrawHitCollection const& strawhits) const { return strawhits[_hind1]; }
//      StrawHit const& sh2(StrawHitCollection const& strawhits) const { return strawhits[_hind2]; }
//      Straw const& s1(StrawHitCollection const& strawhits,Tracker const& tracker) const
//	{ return tracker.getStraw(sh1(strawhits).strawIndex()); }
//      Straw const& s2(StrawHitCollection const& strawhits,Tracker const& tracker) const
//	{ return tracker.getStraw(sh2(strawhits).strawIndex()); }
//      CLHEP::Hep3Vector const& pos() const { return _pos; }
//      // compute the positions for the 2 hits given a direction vector for the track between them.
//      void position(StrawHitCollection const& strawhits,Tracker const& tracker,
//	  CLHEP::Hep3Vector& p1, CLHEP::Hep3Vector& p2, CLHEP::Hep3Vector const& dir) const;
//      float dist() const { return _dist; }
//      PanelId::isep panelSeparation() const { return static_cast<PanelId::isep>(_isep); }
//      float time() const { return _time; }
//      float dt() const { return _dt; } // signed t2 -t1
//      float energy() const { return _edep; }
//      float wdist1() const { return _wd1; }
//      float wdist2() const { return _wd2; }
//      float wdot() const { return _wdot; }
//      std::size_t hitIndex1() const { return _hind1; }
//      std::size_t hitIndex2() const { return _hind2; }
//      float chisq() const { return _chisq; }
//      float mvaout() const { return _mvaout; }
//      void setChisquared(double chisq) { _chisq = chisq; }
//      void setMVAOut(double mvaout) { _mvaout = mvaout; }
//    private:
//      std::size_t _hind1, _hind2; // indices into the straw hit container for the 2 hits making up this stereo hit
//      CLHEP::Hep3Vector _pos; // position in tracker coordinates
//      float _dist; // transverse separation between the 2 hit wires at their POCA
//      int _isep; // separation of measurement planes
//      float _time; // average time of the 2 hits
//      float _dt; // time separation of the 2 hits
//      float _edep; // average energy deposition of the 2 hits
//      float _wd1, _wd2; // distance along the wire of the POCA for the two hits
//      float _chisq; // chisquared of the difference between the stereo position and the time division measurement
//      float _wdot; // dot product of the angle between the 2 wire directions
//      float _mvaout; // output of the MVA
//    };
//   typedef std::vector<mu2e::StereoHit> StereoHitCollection;
//}
//#endif
//
