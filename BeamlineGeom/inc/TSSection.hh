#ifndef BeamlineGeom_TSSection_hh
#define BeamlineGeom_TSSection_hh

//
// Abstract base class for TS sections
//
#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TSSection {

  public:

    virtual ~TSSection(){}

    TSSection() :
      _origin(CLHEP::Hep3Vector()), _rotation(CLHEP::HepRotation()), _materialName("")
    {}

    explicit TSSection(CLHEP::Hep3Vector origin, 
                       CLHEP::HepRotation rotation = CLHEP::HepRotation(),
                       std::string materialName="") :
      _origin(origin),  _rotation(rotation),  _materialName(materialName)
    {}

    CLHEP::Hep3Vector  const & getGlobal()   const { return _origin; }
    CLHEP::HepRotation const * getRotation() const { return &_rotation; }
    std::string const & getMaterial() const { return _materialName; }
    void setMaterial( std::string material ) { _materialName = material; }

  protected:

    CLHEP::Hep3Vector  _origin;
    CLHEP::HepRotation _rotation;
    std::string _materialName;

  };

}
#endif /* BeamlineGeom_TSSection_hh */
