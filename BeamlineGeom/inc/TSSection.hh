#ifndef BeamlineGeom_TSSection_hh
#define BeamlineGeom_TSSection_hh

//
// Abstract base class for TS sections
//
#include <memory>

//#include "BeamlineGeom/inc/TSSection.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TSSection {

  public:

    virtual ~TSSection(){}

    TSSection() :
      _origin(CLHEP::Hep3Vector()), _rotation(CLHEP::HepRotation()) 
    {}

    explicit TSSection(CLHEP::Hep3Vector origin, 
                       CLHEP::HepRotation rotation = CLHEP::HepRotation()) :
      _origin(origin),  _rotation(rotation)
    {}

    virtual CLHEP::Hep3Vector  const & getGlobal()   const { return _origin; }
    virtual CLHEP::HepRotation const * getRotation() const { return &_rotation; }

  protected:

    CLHEP::Hep3Vector  _origin;
    CLHEP::HepRotation _rotation;

  };

}
#endif /* BeamlineGeom_TSSection_hh */
