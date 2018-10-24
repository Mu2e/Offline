////
//// simple reconstruction of 2 straw hits in different views giving 3-d information from stereo 
////
//// $Author: brownd $
//// $Date: 2013/08/21 23:33:41 $
////
//// Original author David Brown
////
//
//// Mu2e includes
//#include "RecoDataProducts/inc/StereoHit.hh"
//#include "Mu2eUtilities/inc/TwoLinePCA.hh"
//// C++ includes
//#include <ostream>
//// clhep includes
//#include "CLHEP/Matrix/Matrix.h"
//#include "CLHEP/Matrix/Vector.h"
//// Framework includes.
//#include "cetlib_except/exception.h"
//
//using namespace std;
//using CLHEP::Hep3Vector;
//using CLHEP::HepVector;
//using CLHEP::HepMatrix;
//
//namespace mu2e {
//  StereoHit::StereoHit() : _hind1(0), _hind2(0), _dist(-5.0), _isep(PanelId::apart), _time(0.0), _dt(0.0), _edep(0.0), _wd1(0.0), _wd2(0.0), _chisq(-1.0)
//  {}
//
//  StereoHit::StereoHit(StrawHitCollection const& strawhits,Tracker const& tracker, size_t ind1, size_t ind2) : 
//    _hind1(ind1),_hind2(ind2), _dist(-5.0), _chisq(-1.0)
//  {
//    // find the midpoint and wire direction of these straws
//    const Hep3Vector& wpos1 = s1(strawhits,tracker).getMidPoint();
//    const Hep3Vector& wdir1 = s1(strawhits,tracker).getDirection();
//    const Hep3Vector& wpos2 = s2(strawhits,tracker).getMidPoint();
//    const Hep3Vector& wdir2 = s2(strawhits,tracker).getDirection();
//    TwoLinePCA pca(wpos1,wdir1,wpos2,wdir2);
//    if(pca.closeToParallel()){  
//      cet::exception("RECO")<<"mu2e::StereoHit: parallel wires" << std::endl;
//    }
//    _pos = 0.5*(pca.point1() + pca.point2());
//    _dist = pca.dca();
//    double t1 = sh1(strawhits).time();
//    double t2 = sh2(strawhits).time();
//    _time = 0.5*(t1+t2);
//    _dt = t2 -t1;
//    _edep = 0.5*(sh1(strawhits).energyDep() + sh2(strawhits).energyDep());
//    _wd1 = pca.s1();
//    _wd2 = pca.s2();
//    _wdot = wdir1.dot(wdir2);
//    _isep = s1(strawhits,tracker).id().separation(s2(strawhits,tracker).id());
//  }
//
//  void StereoHit::position(StrawHitCollection const& strawhits,Tracker const& tracker,
//      Hep3Vector& p1, Hep3Vector& p2, Hep3Vector const& dir) const {
//  // preset to the default position
//      p1 = _pos; p2 = _pos;
//// quick check: if dir == z axis, this is the same as the default calculation
//    if(fabs(dir.z()) == 1.0) return;
//// get straw positions and wire directions
//    const Hep3Vector& wpos1 = s1(strawhits,tracker).getMidPoint();
//    const Hep3Vector& wdir1 = s1(strawhits,tracker).getDirection();
//    const Hep3Vector& wpos2 = s2(strawhits,tracker).getMidPoint();
//    const Hep3Vector& wdir2 = s2(strawhits,tracker).getDirection();
////create linear algebra objects
//    double alpha = (wpos2.z()-wpos1.z())/dir.z();
//    HepVector beta(2);
//    beta(1) = alpha*dir.x() + wpos1.x() - wpos2.x(); 
//    beta(2) = alpha*dir.y() + wpos1.y() - wpos2.y(); 
//    HepMatrix omega(2,2);
//    omega(1,1) = -wdir1.x();
//    omega(1,2) = wdir2.x();
//    omega(2,1) = -wdir1.y();
//    omega(2,2) = wdir2.y();
//    // solve
//    HepVector len = qr_solve(omega,beta);
//    // positions are now given WRT the wire center
//    p1 = wpos1 + len(1)*wdir1;
//    p2 = wpos2 + len(2)*wdir2;
//  }
//
//}
//
//
//
//
