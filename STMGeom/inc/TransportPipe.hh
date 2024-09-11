#ifndef STMGeom_TransportPipe_hh
#define STMGeom_TransportPipe_hh

// Transport Pipe Object
//
// Author: Anthony Palladino
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TransportPipe {
  public:
    TransportPipe(bool build, double radiusIn, double radiusOut,
                  std::string const& material, std::string const& gasMaterial,
                  double upStrSpace, double dnStrHalflength,
                  std::string const& upStrWindowMaterial, double upStrWindowHalflength,
                  std::string const& dnStrWindowMaterial, double dnStrWindowHalflength,
                  double flangeHalflength, double flangeOverhangR,
                  CLHEP::Hep3Vector const& originInMu2e, CLHEP::HepRotation const& rotation
                ) :
      _build( build ),
      _radiusIn( radiusIn ),
      _radiusOut( radiusOut ),
      _material( material ),
      _gasMaterial( gasMaterial ),
      _upStrSpace( upStrSpace ),
      _dnStrHalflength( dnStrHalflength ),
      _upStrWindowMaterial( upStrWindowMaterial ),
      _upStrWindowHalflength( upStrWindowHalflength ),
      _dnStrWindowMaterial( dnStrWindowMaterial ),
      _dnStrWindowHalflength( dnStrWindowHalflength ),
      _flangeHalflength( flangeHalflength ),
      _flangeOverhangR( flangeOverhangR ),
      _originInMu2e( originInMu2e ),
      _rotation    ( rotation     )
    {}

    bool   build()                            const { return _build; }
    double radiusIn()                         const { return _radiusIn; }
    double radiusOut()                        const { return _radiusOut; }
    std::string const & material()            const { return _material; }
    std::string const & gasMaterial()         const { return _gasMaterial; }
    double upStrSpace()                       const { return _upStrSpace; }
    double dnStrHalflength()                  const { return _dnStrHalflength; }
    std::string const & upStrWindowMaterial() const { return _upStrWindowMaterial; }
    double upStrWindowHalflength()            const { return _upStrWindowHalflength; }
    std::string const & dnStrWindowMaterial() const { return _dnStrWindowMaterial; }
    double dnStrWindowHalflength()            const { return _dnStrWindowHalflength; }
    double flangeHalfLength()                 const { return _flangeHalflength; }
    double flangeOverhangR()                  const { return _flangeOverhangR; }

    //double zBegin()          const { return _originInMu2e.z() - zTabletopHalflength(); }
    //double zEnd()            const { return _originInMu2e.z() + zTabletopHalflength(); }

    CLHEP::Hep3Vector const &  originInMu2e() const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()     const { return _rotation; }

    // Genreflex can't do persistency of vector<TransportPipe> without a default constructor
    TransportPipe() {}
  private:

    bool               _build;
    double             _radiusIn;
    double             _radiusOut;
    std::string        _material;
    std::string        _gasMaterial;
    double             _upStrSpace;
    double             _dnStrHalflength;
    std::string        _upStrWindowMaterial;
    double             _upStrWindowHalflength;
    std::string        _dnStrWindowMaterial;
    double             _dnStrWindowHalflength;
    double             _flangeHalflength;
    double             _flangeOverhangR;
    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume

  };

}

#endif/*STMGeom_TransportPipe_hh*/
