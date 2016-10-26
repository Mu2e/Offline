#ifndef STMGeom_ShieldPipe_hh
#define STMGeom_ShieldPipe_hh

// Shield Pipe Object
//
// Author: Anthony Palladino
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class ShieldPipe {
  public:
    ShieldPipe(bool build, double radiusIn, double linerWidth, double radiusOut, double pipeHalfLength,
               std::string materialLiner, std::string material,
               double upStrSpace, double dnStrSpace, 
               double dnStrWallHalflength,
               CLHEP::Hep3Vector originInMu2e, CLHEP::HepRotation rotation
              ) :
      _build( build ),   
      _radiusIn( radiusIn ),
      _linerWidth( linerWidth ),
      _radiusOut( radiusOut ),
      _pipeHalfLength( pipeHalfLength ),
      _materialLiner( materialLiner ),
      _material( material ),
      _upStrSpace( upStrSpace ),
      _dnStrSpace( dnStrSpace ),
      _dnStrWallHalflength( dnStrWallHalflength ),
      _originInMu2e( originInMu2e ),
      _rotation    ( rotation     )
    {}

    bool   build()                            const { return _build; }
    double radiusIn()                         const { return _radiusIn; }
    double linerWidth()                       const { return _linerWidth; }
    double radiusOut()                        const { return _radiusOut; }
    double pipeHalfLength()                   const { return _pipeHalfLength; }
    std::string const & materialLiner()       const { return _materialLiner; }
    std::string const & material()            const { return _material; }      
    double upStrSpace()                       const { return _upStrSpace; }
    double dnStrSpace()                       const { return _dnStrSpace; }
    double dnStrWallHalflength()              const { return _dnStrWallHalflength; }
    CLHEP::Hep3Vector const &  originInMu2e() const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()     const { return _rotation; }
    //double zBegin()          const { return _originInMu2e.z() - zTabletopHalflength(); }
    //double zEnd()            const { return _originInMu2e.z() + zTabletopHalflength(); }

    // Genreflex can't do persistency of vector<ShieldPipe> without a default constructor
    ShieldPipe() {}
  private:

    bool               _build;
    double             _radiusIn;
    double             _linerWidth;
    double             _radiusOut;
    double             _pipeHalfLength;
    std::string        _materialLiner;
    std::string        _material;
    double             _upStrSpace;
    double             _dnStrSpace;
    double             _dnStrWallHalflength;
    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume

  };

}

#endif/*STMGeom_ShieldPipe_hh*/
