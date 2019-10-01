#ifndef GeneralUtilities_HepTransform_HH
#define GeneralUtilities_HepTransform_HH
// 
// A class to represent a rotation about local axes
// followed by a displacement.  Operates on Hep3Vector.
//
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include <ostream>

namespace mu2e {

  class HepTransform {
  public:
    HepTransform():_displacement(0,0,0),
		   _rotation(CLHEP::HepRotation::IDENTITY) {}
    HepTransform( CLHEP::Hep3Vector & disp, CLHEP::HepRotation & rot ) :
      _displacement(disp),
      _rotation(rot) {}
    HepTransform( double dx, double dy, double dz,  
		  double rx, double ry, double rz):
      _displacement(dx,dy,dz) {
      _rotation.rotateX(rx);
      _rotation.rotateY(ry);
      _rotation.rotateZ(rz);
    }

    CLHEP::Hep3Vector  displacement() const { return _displacement;}
    CLHEP::HepRotation rotation()     const { return _rotation;    }

    void  setDisplacement ( CLHEP::Hep3Vector & aDisp ) { 
                                   _displacement = aDisp;}
    void  setRotation ( CLHEP::HepRotation & aRot ) { _rotation = aRot; }

    // combine HepTransform
    HepTransform& operator*=(HepTransform const& b);
    friend HepTransform operator*(HepTransform const& a, HepTransform const& b);

    // apply to a vector
    friend CLHEP::Hep3Vector operator*(HepTransform const& a, 
				CLHEP::Hep3Vector const& v);

  private:
    CLHEP::Hep3Vector  _displacement;
    CLHEP::HepRotation _rotation;
  };

  std::ostream& operator<<(std::ostream& os, const HepTransform& rhs );

} // end of namespace mu2e

#endif  //  GeneralUtilties_HepTransform_HH
