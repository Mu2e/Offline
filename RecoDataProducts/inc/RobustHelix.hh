//
//  Parameters for a purely geometric helix that has non-zero pitch and radius
//  Original Author Dave Brown (LBNL) 15/8/2016
//
#ifndef RecoDataProducts_RobustHelix_HH
#define RecoDataProducts_RobustHelix_HH
// Mu2e
#include "DataProducts/inc/Helicity.hh"
#include "DataProducts/inc/XYZVec.hh"
// Root
#include "Rtypes.h"
// C++
#include <ostream>
namespace mu2e {
  class HelixVal;
  struct RobustHelix {
    RobustHelix(Float_t rcent, Float_t fcent, Float_t radius, Float_t lambda, Float_t fz0) :
      _rcent(rcent), _fcent(fcent), _radius(radius), _lambda(lambda), _fz0(fz0) {}
    RobustHelix() : _rcent(-1.0), _fcent(0.0), _radius(-1.0), _lambda(0.0) , _fz0(0.0), _helicity(Helicity::unknown) {}
    ~RobustHelix(){}
    // accessors
    Float_t radius() const { return _radius; } 
    Float_t rcent() const { return _rcent; } 
    Float_t fcent() const { return _fcent; } 
    Float_t lambda() const { return _lambda; }
    Float_t fz0() const { return _fz0; }
    // simple functions that can be derived from the data members
    // scalar 'momentum' in units of mm
    float momentum() const { return sqrt(_radius*_radius+_lambda*_lambda); }
    // the following is signed by Pz (pitch = dz/dL)
    float pitch() const { return copysign(1.0,_radius)*_lambda/sqrt(_radius*_radius+_lambda*_lambda); }
    Helicity helicity() const { return _helicity; }
    float centerx() const { return _rcent*cos(_fcent); }
    float centery() const { return _rcent*sin(_fcent); }
    XYZVec center() const { return XYZVec(centerx(),centery(),0.0); }
    CLHEP::Hep3Vector centerCLHEP() const { return CLHEP::Hep3Vector(centerx(),centery(),0.0); }
    float chi2dXY  () const { return _chi2dXY; }
    float chi2dZPhi() const { return _chi2dZPhi; }
    // azimuth wrt the circle center expected for a given z position
    Float_t circleAzimuth( double zpos) const { return _lambda != 0.0 ? _fz0 + zpos/_lambda : 0.0; }
    // position in space given the Z position of the input vector
    void position(XYZVec& pos) const {
      pos.SetX(centerx() + _radius*cos(circleAzimuth(pos.z())));
      pos.SetY(centery() + _radius*sin(circleAzimuth(pos.z())));
    //  pos.setz(0.0);  not sure why this was here
    }
    // unit vector in direction at the given z
    void direction(float zval,XYZVec& dir) const {
      float mom = momentum();
      dir.SetX(-_radius*sin(circleAzimuth(zval))/mom);
      dir.SetY(_radius*cos(circleAzimuth(zval))/mom);
      dir.SetZ(_lambda/mom);
    }
    void position(CLHEP::Hep3Vector& pos) const {
      pos.setX(centerx() + _radius*cos(circleAzimuth(pos.z())));
      pos.setY(centery() + _radius*sin(circleAzimuth(pos.z())));
    //  pos.setz(0.0);  not sure why this was here
    }
    // unit vector in direction at the given z
    void direction(float zval,CLHEP::Hep3Vector& dir) const {
      float mom = momentum();
      dir.setX(-_radius*sin(circleAzimuth(zval))/mom);
      dir.setY(_radius*cos(circleAzimuth(zval))/mom);
      dir.setZ(_lambda/mom);
    }
    // test self-consistency
    bool validHelicity() const { return Helicity(_lambda) == _helicity; }
    // data members
    Float_t _rcent; // radius of circle center
    Float_t _fcent; // azimuth of circle center
    Float_t _radius;  // transverse radius of the helix (mm).  Always positive
    Float_t _lambda; // dz/dphi (mm/radian)
    Float_t _fz0; // azimuth (phi) at the center z position (radians)
    Float_t _chi2dXY;   // squared sum of the hit residuals wrt the XY plane normalized by the degree of freedom
    Float_t _chi2dZPhi; // squared sum of the hit residuals wrt the ZPhi plane normalized by the degree of freedom
    Helicity _helicity; // explicit helicity value
  };
  std::ostream& operator<<(std::ostream& os, mu2e::RobustHelix const& helix);

}
#endif
