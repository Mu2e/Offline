#ifndef BeamlineGeom_StraightSection_hh
#define BeamlineGeom_StraightSection_hh

//
// Class to represent the transport solenoid
//
#include <memory>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class StraightSection {

  friend class BeamlineMaker;

  public:

    StraightSection() { _rotation=0;};

    explicit StraightSection(double halfZ, CLHEP::Hep3Vector origin, CLHEP::HepRotation *rotation) :
      _halfZ(halfZ), _origin(origin), _rotation(rotation) {;};

    ~StraightSection(){ if( _rotation!=0 ) delete _rotation; };

    void set(double halfZ, CLHEP::Hep3Vector origin, CLHEP::HepRotation *rotation) {
      _halfZ=halfZ;
      _origin=origin;
      _rotation=rotation;
    };

    double getHalfLength() const { return _halfZ; }
    CLHEP::Hep3Vector const& getGlobal() const { return _origin; }
    CLHEP::HepRotation * getRotation() const { return _rotation; }

  private:

    double _halfZ;
    CLHEP::Hep3Vector _origin;
    CLHEP::HepRotation * _rotation;
};

}
#endif /* BeamlineGeom_StraightSection_hh */
