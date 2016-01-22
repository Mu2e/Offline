#ifndef STMGeom_PermanentMagnet_hh
#define STMGeom_PermanentMagnet_hh

// Permanent Magnet Object
//
// Author: Anthony Palladino
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class PermanentMagnet {
  public:
    PermanentMagnet(bool build,
                    double halfX, double halfY, double halfZ,
                    double holeHalfX, double holeHalfY, 
                    CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(),
                    CLHEP::HepRotation const & rotation = CLHEP::HepRotation(), 
                    std::string const & materialName = "", double field=0.0,
                    bool fieldVisible=false
                   ) :
      _build( build),
      _xhl( halfX ),
      _yhl( halfY ),
      _zhl( halfZ ),
      _holexhl( holeHalfX ),
      _holeyhl( holeHalfY ),
      _originInMu2e( originInMu2e ),
      _rotation    ( rotation     ),
      _materialName( materialName ),
      _field( field ),
      _fieldVisible ( fieldVisible )
    {}

    bool   build()           const { return _build; }
    double xHalfLength()     const { return _xhl; }
    double yHalfLength()     const { return _yhl; }
    double zHalfLength()     const { return _zhl; }
    double xHoleHalfLength() const { return _holexhl; }
    double yHoleHalfLength() const { return _holeyhl; }
    double zBegin()          const { return _originInMu2e.z() - zHalfLength(); }
    double zEnd()            const { return _originInMu2e.z() + zHalfLength(); }
    double field()           const { return _field; }
    bool   fieldVisible()    const { return _fieldVisible; }
   
    CLHEP::Hep3Vector const &  originInMu2e() const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()     const { return _rotation; }
    std::string const &        materialName() const { return _materialName; }    

    // Genreflex can't do persistency of vector<PermanentMagnet> without a default constructor
    PermanentMagnet() {}
  private:

    bool   _build;
    double _xhl;
    double _yhl;
    double _zhl;
    double _holexhl;
    double _holeyhl;
    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume
    std::string        _materialName;
    double _field;
    bool   _fieldVisible;

  };

}

#endif/*STMGeom_PermanentMagnet_hh*/
