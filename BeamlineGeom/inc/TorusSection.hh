#ifndef BeamlineGeom_TorusSection_hh
#define BeamlineGeom_TorusSection_hh

//
// Class to represent the transport solenoid
//
#include <memory>

#include "BeamlineGeom/inc/TSSection.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TorusSection : public TSSection {

  friend class BeamlineMaker;

  public:

    TorusSection() : 
      _rTorus(0.),_rIn(0.),_rOut(0.),
      _phiBegin(0.),_deltaPhi(0.) 
    {
      _rotation.reset( nullptr );  
    }

    TorusSection(double rTorus, double rIn, double rOut, double phi0, double dPhi, 
                 CLHEP::Hep3Vector origin, CLHEP::HepRotation *rotation) :
      _rTorus(rTorus),_rIn(rIn),_rOut(rOut),
      _phiBegin(phi0),_deltaPhi(dPhi)
    {
      _origin=origin; // _origin is derived data member; cannot be in initialization list
      _rotation.reset( rotation ) ;
    }

    ~TorusSection(){}

    void set(double rTorus, double rIn, double rOut, double phi0, double dPhi, 
             CLHEP::Hep3Vector origin, CLHEP::HepRotation *rotation) {
      _rTorus=rTorus;
      _rIn=rIn;
      _rOut=rOut;
      _phiBegin=phi0;
      _deltaPhi=dPhi;
      _origin=origin;
      _rotation.reset( rotation );
    }

    double torusRadius() const {return _rTorus;}
    double rIn()         const {return _rIn;   }
    double rOut()        const {return _rOut;  }
    double phiStart()    const {return _phiBegin; }
    double deltaPhi()    const {return _deltaPhi; }

  private:

    // All dimensions in mm.
    double _rTorus;
    double _rIn;
    double _rOut;

    double _phiBegin;
    double _deltaPhi;

};

}
#endif /* BeamlineGeom_TorusSection_hh */
