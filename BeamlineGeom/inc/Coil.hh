#ifndef BeamlineGeom_Coil_hh
#define BeamlineGeom_Coil_hh

//
// Class to represent the transport solenoid
//

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class Coil {

  public:

    Coil() : _rIn(0.), _rOut(0.), _halfZ(0.), _origin(), _rotation() {}

    Coil( double x, double y, double z,
          double rIn, double rOut, double halfZ,
          CLHEP::HepRotation const& rotation = CLHEP::HepRotation() ) :
      _rIn( rIn ), _rOut( rOut ) , _halfZ( halfZ ),
      _origin( x, y, z ),
      _rotation(rotation)
    {
    }

    double rIn()        const { return _rIn;   }
    double rOut()       const { return _rOut;  }
    double halfLength() const { return _halfZ; }
    CLHEP::Hep3Vector  const & getGlobal()   const { return _origin; }
    CLHEP::HepRotation const * getRotation() const { return &_rotation; }

    void  setRotation( CLHEP::HepRotation const& rotation ) { _rotation = rotation; }

    void setOrigin  ( const double x, const double y, const double z ) {
      _origin = CLHEP::Hep3Vector( x, y, z );
    }

  private:

    double _rIn;
    double _rOut;
    double _halfZ;
    CLHEP::Hep3Vector _origin;
    CLHEP::HepRotation _rotation;

  };

}
#endif /* BeamlineGeom_Coil_hh */
