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

  friend class BeamlineMaker;

  public:

    Coil() : _rIn(0.), _rOut(0.), _halfZ(0.) {}
    
    Coil( double rIn, double rOut, double halfZ, CLHEP::HepRotation * rotation = nullptr ) :
      _rIn( rIn ), _rOut( rOut ) , _halfZ( halfZ ) 
    {
      _rotation.reset( rotation );
    }

    double rIn()        const { return _rIn;   }
    double rOut()       const { return _rOut;  }
    double halfLength() const { return _halfZ; }
    CLHEP::Hep3Vector const& getGlobal() const { return _origin; }
    CLHEP::HepRotation * getRotation() const { return _rotation.get(); }

    void  setRotation( CLHEP::HepRotation * rotation ) { _rotation.reset( rotation ); }

    inline void setOrigin  ( const unsigned nCoils, 
                             TSSection*    ts,
                             const unsigned iCoil,
                             const CLHEP::HepRotation * rotMatrix,
                             const double totCoilL ) {

      TorusSection* torsec    ( dynamic_cast<TorusSection*>( ts ) );
      StraightSection* strsec ( dynamic_cast<StraightSection*>( ts ) ); 

      static double str_offset{};

      if ( torsec ) {
        str_offset = 0.;
        // don't use total coil length here; spacing determined by angle
        const double R           = torsec->torusRadius();
        const double thetaOffset = atan( this->halfLength()/R );
        const double thetaIncr   = (CLHEP::halfpi - 2*thetaOffset)/(nCoils - 1);
        const double sgn    = std::copysign(1.0,torsec->getGlobal().x() ); 
        const double theta  = sgn > 0 ? // check if torus is in upstream/downstream portion 
          -(thetaOffset + iCoil*thetaIncr ) :
          3*CLHEP::halfpi-(thetaOffset + iCoil*thetaIncr );

        CLHEP::Hep3Vector pos( R*cos(theta), torsec->getGlobal().y(), -sgn*R*sin(theta) );
        _origin = CLHEP::Hep3Vector( pos + torsec->getGlobal() );
      }
      if ( strsec ) {
        const double spacehl = (strsec->getHalfLength() - totCoilL*.5)/(nCoils+1);
        if ( iCoil == 0 ) str_offset -= strsec->getHalfLength();
        str_offset += 2*spacehl + this->halfLength();
        CLHEP::Hep3Vector pos ( 0., 0., str_offset );
        _origin = CLHEP::Hep3Vector( (*rotMatrix)*pos + strsec->getGlobal() );
        str_offset += this->halfLength();
      }
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
