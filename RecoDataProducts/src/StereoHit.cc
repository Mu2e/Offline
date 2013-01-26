//
// simple reconstruction of 2 straw hits in different views giving 3-d information from stereo 
//
// $Author: brownd $
// $Date: 2013/01/26 18:18:44 $
//
// Original author David Brown
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib/exception.h"

// Mu2e includes
#include "RecoDataProducts/inc/StereoHit.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
using namespace std;

namespace mu2e {
  StereoHit::StereoHit() : _hind1(0), _hind2(0), _dist(0.0), _time(0.0), _dt(0.0), _wd1(0.0), _wd2(0.0)
  {}

  StereoHit::StereoHit(StrawHitCollection const& strawhits,Tracker const& tracker, size_t ind1, size_t ind2) : 
    _hind1(ind1),_hind2(ind2), _dist(-1.0)
  {
    // find the midpoint and wire direction of these straws
    const CLHEP::Hep3Vector& wpos1 = s1(strawhits,tracker).getMidPoint();
    const CLHEP::Hep3Vector& wdir1 = s1(strawhits,tracker).getDirection();
    const CLHEP::Hep3Vector& wpos2 = s2(strawhits,tracker).getMidPoint();
    const CLHEP::Hep3Vector& wdir2 = s2(strawhits,tracker).getDirection();
    TwoLinePCA pca(wpos1,wdir1,wpos2,wdir2);
    if(pca.closeToParallel()){  
      cet::exception("RECO")<<"mu2e::StereoHit: parallel wires" << std::endl;
    }
    _pos = 0.5*(pca.point1() + pca.point2());
    _dist = pca.dca();
    double t1 = sh1(strawhits).time();
    double t2 = sh2(strawhits).time();
    _time = 0.5*(t1+t2);
    _dt = t2 -t1;
    _wd1 = pca.s1();
    _wd2 = pca.s2();
  }
}


