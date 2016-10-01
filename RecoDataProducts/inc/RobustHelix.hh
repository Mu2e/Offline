//
//  Parameters for a purely geometric helix that has non-zero pitch and radius
//  Original Author Dave Brown (LBNL) 15/8/2016
//
#ifndef RecoDataProducts_RobustHelix_HH
#define RecoDataProducts_RobustHelix_HH
// Mu2e
#include "DataProducts/inc/Helicity.hh"
// CLHEP
#include "CLHEP/Vector/ThreeVector.h"
// Root
#include "Rtypes.h"
// C++
#include <ostream>
namespace mu2e {
  class HelixVal;
  class RobustHelix {
    public:
      RobustHelix(CLHEP::Hep3Vector const& center, Float_t radius, Float_t lambda, Float_t fz0) :
	_center(center), _radius(radius), _lambda(lambda), _fz0(fz0) {}
      RobustHelix() : _radius(-1.0), _lambda(0.0) , _fz0(0.0) {}
      ~RobustHelix(){}
      // accessors
      CLHEP::Hep3Vector const& center() const { return _center; }
      Float_t radius() const { return _radius; } 
      Float_t lambda() const { return _lambda; }
      Float_t fz0() const { return _fz0; }
      Helicity helicity() const { return Helicity(static_cast<float>(_lambda)); }
      // non-const accessors for setting values
      CLHEP::Hep3Vector& center() { return _center; }
      Float_t& radius() { return _radius; } 
      Float_t& lambda() { return _lambda; }
      Float_t& fz0() { return _fz0; }
      // azimuth wrt the circle center expected for a given z position
      Float_t circleAzimuth( double zpos) const { return _fz0 + (zpos - _center.z())/_lambda; }
      // position in space given the Z position of the input vector
      void position(CLHEP::Hep3Vector& pos) const {
	pos.setX(_center.x() + _radius*cos(circleAzimuth(pos.z())));
	pos.setY(_center.y() + _radius*sin(circleAzimuth(pos.z())));
      }
    private:
      CLHEP::Hep3Vector _center; // circle center.  The z coordinate is not a fit parameter (nominally 0) (mm)
      Float_t _radius;  // transverse radius of the helix (mm).  Always positive
      Float_t _lambda; // dz/dphi (mm/radian)
      Float_t _fz0; // azimuth (phi) at the center z position (radians)
  };
  std::ostream& operator<<(std::ostream& os, mu2e::RobustHelix const& helix);

}
#endif
