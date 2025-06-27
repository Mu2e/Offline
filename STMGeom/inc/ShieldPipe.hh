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
    ShieldPipe(bool build, double radiusIn, bool hasLiner, double linerWidth, double radiusOut, double pipeHalfLength,
               std::string const& materialLiner, std::string const& material,
               bool matchPipeBlock,
               double upStrSpace, double dnStrSpace,
               double dnStrWallHalflength, double dnStrWallHoleRadius,
               double dnStrWallHalfHeight, double dnStrWallHalfWidth,
               double dnStrWallGap, std::string const& dnStrWallMaterial, bool buildMatingBlock,
               double upStrAirGap,
               CLHEP::Hep3Vector const& originInMu2e, CLHEP::HepRotation const& rotation
              ) :
      _build( build ),
      _radiusIn( radiusIn ),
      _hasLiner( hasLiner ),
      _linerWidth( linerWidth ),
      _radiusOut( radiusOut ),
      _pipeHalfLength( pipeHalfLength ),
      _materialLiner( materialLiner ),
      _material( material ),
      _matchPipeBlock( matchPipeBlock ),
      _upStrSpace( upStrSpace ),
      _dnStrSpace( dnStrSpace ),
      _dnStrWallHalflength( dnStrWallHalflength ),
      _dnStrWallHoleRadius( dnStrWallHoleRadius ),
      _dnStrWallHalfHeight( dnStrWallHalfHeight ),
      _dnStrWallHalfWidth( dnStrWallHalfWidth ),
      _dnStrWallGap( dnStrWallGap ),
      _dnStrWallMaterial( dnStrWallMaterial ),
      _buildMatingBlock( buildMatingBlock ),
      _upStrAirGap( upStrAirGap ),
      _originInMu2e( originInMu2e ),
      _rotation    ( rotation     )
    {}

    bool   build()                            const { return _build; }
    bool   hasLiner()                         const { return _hasLiner; }
    double radiusIn()                         const { return _radiusIn; }
    double linerWidth()                       const { return _linerWidth; }
    double radiusOut()                        const { return _radiusOut; }
    double pipeHalfLength()                   const { return _pipeHalfLength; }
    std::string const & materialLiner()       const { return _materialLiner; }
    std::string const & material()            const { return _material; }
    bool   matchPipeBlock()                   const { return _matchPipeBlock; }
    double upStrSpace()                       const { return _upStrSpace; }
    double dnStrSpace()                       const { return _dnStrSpace; }
    double dnStrWallHalflength()              const { return _dnStrWallHalflength; }
    double dnStrWallHoleRadius()              const { return _dnStrWallHoleRadius; }
    double dnStrWallHalfHeight()              const { return _dnStrWallHalfHeight; }
    double dnStrWallHalfWidth()               const { return _dnStrWallHalfWidth; }
    double dnStrWallGap()                     const { return _dnStrWallGap; }
    std::string const & dnStrWallMaterial()   const { return _dnStrWallMaterial; }
    bool buildMatingBlock()                   const { return _buildMatingBlock; }
    double upStrAirGap()                       const { return _upStrAirGap; }
    CLHEP::Hep3Vector const &  originInMu2e() const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()     const { return _rotation; }
    //double zBegin()          const { return _originInMu2e.z() - zTabletopHalflength(); }
    //double zEnd()            const { return _originInMu2e.z() + zTabletopHalflength(); }

    // Genreflex can't do persistency of vector<ShieldPipe> without a default constructor
    ShieldPipe() {}
  private:

    bool               _build;
    double             _radiusIn;
    bool               _hasLiner;
    double             _linerWidth;
    double             _radiusOut;
    double             _pipeHalfLength;
    std::string        _materialLiner;
    std::string        _material;
    bool               _matchPipeBlock;
    double             _upStrSpace;
    double             _dnStrSpace;
    double             _dnStrWallHalflength;
    double             _dnStrWallHoleRadius;
    double             _dnStrWallHalfHeight;
    double             _dnStrWallHalfWidth;
    double             _dnStrWallGap; //between mating block/shield pipe and magnet
    std::string        _dnStrWallMaterial;
    bool               _buildMatingBlock;
    double             _upStrAirGap;
    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume

  };

}

#endif/*STMGeom_ShieldPipe_hh*/
