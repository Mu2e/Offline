#ifndef CIRCLE2D_HH
#define CIRCLE2D_HH

#include "CDFTrajectory/inc/Trajectory.hh"
namespace mu2e {
class Circle2D : public Trajectory  {

public:

  // Constructors
  Circle2D() {}
  inline Circle2D(double a, double phi0, double b);

  // Access parameters
  double getA(void)         const;
  double getB(void)         const;
  double getPhi0(void)      const;
  double getCurvature(void) const;
  double getD0(void)        const;

  // Points on the curve
  double   getPhiAtR(double rho) const;

  // Memory
  Circle2D(const Circle2D &right);
  const Circle2D & operator=(const Circle2D &right);
  ~Circle2D();

  // support for generic fitting
  void setParameters(CLHEP::HepVector pars);
  virtual HepGeom::Point3D<double> getPosition(double s = 0.0) const;
  virtual HepGeom::Vector3D<double> getDirection(double s = 0.0) const;
  virtual HepGeom::Vector3D<double> getSecondDerivative(double s = 0.0) const;
  static Circle2D create(const CLHEP::HepVector & v);
  static unsigned int getParameterSpaceSize();
  CLHEP::HepVector getParameters() const;


private:

  double _a;
  double _b;
  double _phi0;

};


// _a = c/(1.0 + 2.0*_c*_delta);
// _b = _delta - _a*_delta*_delta;
// _delta = 2.0*_b/(1.0 + sqrt(1.0-4.0*_a*_b));
// _c     = _a/(1.0-2.0*_a*_delta);


//Constructors

inline
Circle2D::Circle2D(double a,double phi0,double b) : _a(a),_phi0(phi0),_b(b) {}


// Access Parameters
inline
double Circle2D::getA(void) const {return _a;}

inline
double Circle2D::getB(void) const {return _b;}

inline
double Circle2D::getPhi0(void) const {return _phi0;}

inline
Circle2D::Circle2D(const Circle2D &right)
  : _phi0(right._phi0), _a(right._a), _b(right._b)
{}


// Points on the curve
inline
double Circle2D::getPhiAtR(double rho) const {
  double sinphi = _a*rho + _b/rho;
  double phi = asin(sinphi);
  phi += _phi0;
  return phi;
}


// Memory
inline
const Circle2D & Circle2D::operator=(const Circle2D &right) {
  _phi0=right._phi0;
  _a = right._a;
  _b = right._b;

  return *this;
}

inline
Circle2D::~Circle2D() {}
} // end namespace mu2e

#endif
