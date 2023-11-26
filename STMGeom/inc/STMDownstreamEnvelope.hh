#ifndef STMGeom_STMDownstreamEnvelope_hh
#define STMGeom_STMDownstreamEnvelope_hh

// Mother volume for the STM downstream infrastructure
//
// Author: Andy Edmonds (copued from PermanentMagnet.hh)
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class STMDownstreamEnvelope {
  public:
    STMDownstreamEnvelope(bool build,
                    double halfX, double halfY, double halfZ,
                    CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(),
                    CLHEP::HepRotation const & rotation = CLHEP::HepRotation(),
                    std::string const & materialName = ""
                   ) :
      _build( build),
      _xhl( halfX ),
      _yhl( halfY ),
      _zhl( halfZ ),
      _originInMu2e( originInMu2e ),
      _rotation    ( rotation     ),
      _materialName( materialName )
    {}

    bool   build()           const { return _build; }
    double xHalfLength()     const { return _xhl; }
    double yHalfLength()     const { return _yhl; }
    double zHalfLength()     const { return _zhl; }
    double zBegin()          const { return _originInMu2e.z() - zHalfLength(); }
    double zEnd()            const { return _originInMu2e.z() + zHalfLength(); }

    CLHEP::Hep3Vector const &  originInMu2e() const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()     const { return _rotation; }
    std::string const &        materialName() const { return _materialName; }

    // Genreflex can't do persistency of vector<STMDownstreamEnvelope> without a default constructor
    STMDownstreamEnvelope() {}
  private:

    bool   _build;
    double _xhl;
    double _yhl;
    double _zhl;
    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume
    std::string        _materialName;
  };

}

#endif/*STMGeom_STMDownstreamEnvelope_hh*/
