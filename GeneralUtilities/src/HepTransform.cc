
#include "GeneralUtilities/inc/HepTransform.hh"

using namespace std;

namespace mu2e {

  ostream& operator<<(ostream& os, const HepTransform& rhs) {
    os << "translate: " << rhs.displacement() << endl 
       << "rotate: " 
       << rhs.rotation() << endl;
    return os;
  } // end of outputter


  HepTransform& HepTransform::operator*=(HepTransform const& b) {
    _displacement = _displacement + _rotation * b.displacement();
    _rotation = _rotation * b.rotation();
    return *this;
  }

  HepTransform operator*(HepTransform const& a, HepTransform const& b) {
    CLHEP::Hep3Vector v = a.displacement() + a.rotation() * b.displacement();
    CLHEP::HepRotation r = a.rotation() * b.rotation();
    return HepTransform(v,r);
  }

  CLHEP::Hep3Vector operator*(HepTransform const& a, 
			      CLHEP::Hep3Vector const& v) {
    return a.displacement() + a.rotation()*v;
  }


} // end of namespace mu2e
