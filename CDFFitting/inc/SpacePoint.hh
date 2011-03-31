// The class SpacePoint is basically just a HepPoint3D, dressed
// with a few methods to make it interact with the fitting classes.
//
// Joe Boudreau Sep. 1997.
//
#ifndef SPACEPOINT_HH_
#define SPACEPOINT_HH_
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Matrix/Vector.h"
class SpacePoint: public HepPoint3D {

public:

  // Constructor
  SpacePoint(double x=0, double y=0, double z=0)
    : HepPoint3D(x, y, z) {}
  
  // Copy constructor
  SpacePoint(const Hep3Vector &v) : HepPoint3D(v) {}
  
  // Assignment
  SpacePoint& operator=(const HepPoint3D &v) {
    setX(v.x()); setY(v.y()); setZ(v.z()); return *this;
  }
  
  // Assignment of Hep3Vector and classes derived from it (Vector3D, Normal3D)
  SpacePoint& operator=(const Hep3Vector &v) {
    setX(v.x()); setY(v.y()); setZ(v.z()); return *this;
  }

  bool operator == ( const SpacePoint & right) const {
    return right.x() == x() && right.y() == y() && right.z() == z();
  }

  HepVector getParameters() const {
    HepVector myVector(3);
    myVector(1)=x();
    myVector(2)=y();
    myVector(3)=z();
    return myVector;
  }

  static SpacePoint create(const HepVector & v)  {
    return SpacePoint(v(1),v(2),v(3));
  }

  static unsigned int getParameterSpaceSize() {
    return 3;
  }

private:

};
#endif //SPACEPOINT_HH_
