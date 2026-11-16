#ifndef STMGeom_SSCFrontShield_hh
#define STMGeom_SSCFrontShield_hh

// The shielding in front of the STM spot-size collimator: a wall of
// lead bricks, an aluminium shelf, and two polyethylene blocks, one
// of them bored for the beam.
//
// The bricks themselves are described by LeadBrick, which is shared
// with every other structure that stacks them; this class says only
// which type goes where, in the same BrickWallBrick form the walls
// use. The geometry file groups them by type because a group shares an
// offset rule, but by the time they arrive here they are one list:
// which group a brick came from no longer matters once its offset is
// resolved.
//
// Positions arrive here already resolved. The geometry file writes
// every center as if the structure sat exactly on the SSC axis, and
// STMMaker adds the offsets that displace it -- the whole-structure
// offsetX/offsetY, the LaBr side brick offset, and the poly2 offset --
// so nothing downstream has to know which offset applies to which
// piece.
//
// Author: Yongyi Wu

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/STMGeom/inc/BrickWall.hh"

namespace mu2e {

  class SSCFrontShield {
  public:

    SSCFrontShield(bool build,
                   std::vector<BrickWallBrick> const & bricks,
                   std::string const & shelfMaterial,
                   CLHEP::Hep3Vector const & shelfDim,
                   CLHEP::Hep3Vector const & shelfCenter,
                   std::string const & poly1Material,
                   CLHEP::Hep3Vector const & poly1Dim,
                   CLHEP::Hep3Vector const & poly1Center,
                   double poly1BoreR, double poly1BoreDX, double poly1BoreDY,
                   std::string const & poly2Material,
                   CLHEP::Hep3Vector const & poly2Dim,
                   CLHEP::Hep3Vector const & poly2Center
                   ) :
      _build(build),
      _bricks(bricks),
      _shelfMaterial(shelfMaterial),
      _shelfDim(shelfDim),
      _shelfCenter(shelfCenter),
      _poly1Material(poly1Material),
      _poly1Dim(poly1Dim),
      _poly1Center(poly1Center),
      _poly1BoreR(poly1BoreR),
      _poly1BoreDX(poly1BoreDX),
      _poly1BoreDY(poly1BoreDY),
      _poly2Material(poly2Material),
      _poly2Dim(poly2Dim),
      _poly2Center(poly2Center)
    {
    }

    bool build() const {return _build;}

    // Which brick goes where. Dimensions and material come from
    // LeadBrick. These carry no bores: the beam opening here is a
    // course left empty, not a hole through a brick.
    std::vector<BrickWallBrick> const & bricks() const {return _bricks;}

    // Aluminium shelf, carrying the bricks.
    std::string const & shelfMaterial() const {return _shelfMaterial;}
    CLHEP::Hep3Vector const & shelfDim()    const {return _shelfDim;}
    CLHEP::Hep3Vector const & shelfCenter() const {return _shelfCenter;}

    // The large bored poly block. The bore runs along z through the
    // full depth, offset from the block center by (boreDX, boreDY) so
    // that it lands on the SSC axis.
    std::string const & poly1Material() const {return _poly1Material;}
    CLHEP::Hep3Vector const & poly1Dim()    const {return _poly1Dim;}
    CLHEP::Hep3Vector const & poly1Center() const {return _poly1Center;}
    double poly1BoreR()  const {return _poly1BoreR;}
    double poly1BoreDX() const {return _poly1BoreDX;}
    double poly1BoreDY() const {return _poly1BoreDY;}

    // The small poly block in the beam opening.
    std::string const & poly2Material() const {return _poly2Material;}
    CLHEP::Hep3Vector const & poly2Dim()    const {return _poly2Dim;}
    CLHEP::Hep3Vector const & poly2Center() const {return _poly2Center;}

    // Genreflex can't do persistency of vector<SSCFrontShield> without
    // a default constructor
    SSCFrontShield() {}

  private:

    bool _build;

    std::vector<BrickWallBrick> _bricks;

    std::string _shelfMaterial;
    CLHEP::Hep3Vector _shelfDim;
    CLHEP::Hep3Vector _shelfCenter;

    std::string _poly1Material;
    CLHEP::Hep3Vector _poly1Dim;
    CLHEP::Hep3Vector _poly1Center;
    double      _poly1BoreR;
    double      _poly1BoreDX;
    double      _poly1BoreDY;

    std::string _poly2Material;
    CLHEP::Hep3Vector _poly2Dim;
    CLHEP::Hep3Vector _poly2Center;
  };

}

#endif/*STMGeom_SSCFrontShield_hh*/
