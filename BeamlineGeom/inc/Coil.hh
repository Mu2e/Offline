#ifndef BeamlineGeom_Coil_hh
#define BeamlineGeom_Coil_hh

//
// Class to represent the transport solenoid
//
#include <memory>

#include "BeamlineGeom/inc/TSSection.hh"
#include "BeamlineGeom/inc/TorusSection.hh"
#include "BeamlineGeom/inc/StraightSection.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class Coil {

  public:

    Coil() : _rIn(0.), _rOut(0.), _halfZ(0.) {}
    
    Coil( double x, double y, double z, 
	  double rIn, double rOut, double halfZ, 
	  CLHEP::HepRotation * rotation = nullptr ) :
      _rIn( rIn ), _rOut( rOut ) , _halfZ( halfZ ) 
    {
      _origin = CLHEP::Hep3Vector( x, y, z );
      _rotation.reset( rotation );
    }

    double rIn()        const { return _rIn;   }
    double rOut()       const { return _rOut;  }
    double halfLength() const { return _halfZ; }
    CLHEP::Hep3Vector const& getGlobal() const { return _origin; }
    CLHEP::HepRotation * getRotation() const { return _rotation.get(); }

    void  setRotation( CLHEP::HepRotation * rotation ) { _rotation.reset( rotation ); }

    void setOrigin  ( const double x, const double y, const double z ) {
      _origin = CLHEP::Hep3Vector( x, y, z );
    }

  private:

    double _rIn;
    double _rOut;
    double _halfZ;
    CLHEP::Hep3Vector _origin;
    std::unique_ptr<CLHEP::HepRotation> _rotation;

  };

}
#endif /* BeamlineGeom_Coil_hh */
