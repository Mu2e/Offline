#include "CDFTrajectory/inc/Circle2D.hh"
namespace mu2e {

double Circle2D::getCurvature(void) const { 
  double delta = 2.0*_b/(1.0 + sqrt(1.0-4.0*_a*_b));
  double c     = _a/(1.0-2.0*_a*delta);
  return c;
}

double Circle2D::getD0(void) const { 
  return 2.0*_b/(1.0 + sqrt(1.0-4.0*_a*_b));
}

void Circle2D::setParameters(CLHEP::HepVector pars) {
  _a    = pars[0];
  _b    = pars[1];
  _phi0 = pars[2];

//  _delta = 2.0*_b/(1.0 + sqrt(1.0-4.0*_a*_b));
//  _c     = _a/(1.0-2.0*_a*_delta);
}

  HepGeom::Point3D<double> Circle2D::getPosition(double s ) const {
  std::cerr << "Circle2D not parametrized as a function of arc length" << std::endl;
  HepGeom::Point3D<double> p;
  return p;
}

  HepGeom::Vector3D<double> Circle2D::getDirection(double s ) const {
  std::cerr << "Circle2D not parametrized as a function of arc length" << std::endl;
  HepGeom::Vector3D<double> p;
  return p;
}

  HepGeom::Vector3D<double> Circle2D::getSecondDerivative(double s ) const {
  std::cerr << "Circle2D not parametrized as a function of arc length" << std::endl;
  HepGeom::Vector3D<double> p;
  return p;
}

Circle2D Circle2D::create(const CLHEP::HepVector& v) {return Circle2D(v(1),v(2),v(3));}

unsigned int Circle2D::getParameterSpaceSize() { return 3; }

CLHEP::HepVector Circle2D::getParameters() const {
  CLHEP::HepVector pars(3);
  pars[0] = getA();
  pars[1] = getB();
  pars[2] = getPhi0();
  return pars;
}
} // end namespace mu2e

