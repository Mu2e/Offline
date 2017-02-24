//
// simple reconstruction of 2 straw hits in different views giving 3-d information from stereo 
//
// $Author: brownd $
// $Date: 2013/08/21 23:33:41 $
//
// Original author David Brown
//

// Mu2e includes
#include "RecoDataProducts/inc/StereoHit.hh"
#include "GeneralUtilities/inc/TwoLinePCA.hh"
// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib_except/exception.h"

using namespace std;

namespace mu2e {
  StereoHit::StereoHit() : _hind1(0), _hind2(0), _dist(-5.0), _isep(PanelId::apart), _time(0.0), _dt(0.0), _edep(0.0), _wd1(0.0), _wd2(0.0), _chisq(-1.0)
  {}

  StereoHit::StereoHit(StrawHitCollection const& strawhits,Tracker const& tracker, size_t ind1, size_t ind2) : 
    _hind1(ind1),_hind2(ind2), _dist(-5.0), _chisq(-1.0)
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
    _edep = 0.5*(sh1(strawhits).energyDep() + sh2(strawhits).energyDep());
    _wd1 = pca.s1();
    _wd2 = pca.s2();
    _wdot = wdir1.dot(wdir2);
    _isep = s1(strawhits,tracker).id().getPanelId().separation(s2(strawhits,tracker).id().getPanelId());
  }

}


