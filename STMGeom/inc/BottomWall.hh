#ifndef STMGeom_BottomWall_hh
#define STMGeom_BottomWall_hh

// The floor of the STM shield house, including the steel baseplate the
// whole house stands on.
//
// Five layers worked upward in +y above the plate: a borated poly
// prism, a lead layer, the same prism again, a second lead layer, and
// two copper prisms. 46 pieces, none of them bored.
//
// Both lead layers are mostly courses, which a BrickWall describes, but
// each also carries a few bricks that break the pattern. Those arrive
// already placed in the same BrickWall, since a BrickWallBrick holds
// its own center and orientation -- course brick or stray matters only
// in the config, where one is written as a list and the other as a
// tuple.
//
// Every position arrives resolved from STMMaker. The wall is fixed by
// one reference -- 5 in inboard of the house reference in x, at that
// reference's z, and at the TOP of the baseplate in y.
//
// Author: Yongyi Wu

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/STMGeom/inc/BrickWall.hh"

namespace mu2e {

  // The baseplate. Its own piece rather than a sheet in the layer list,
  // because it is the datum the rest of the house measures from and it
  // sits below the reference while every layer sits above it.
  struct BottomWallPlate {
    std::string       material;
    CLHEP::Hep3Vector halfDim;
    CLHEP::Hep3Vector center;
  };

  // A prism swept along y: the poly L, used twice, and the two copper
  // pieces. Given in the ExtShieldDownstream form, a u/v polygon plus a
  // sweep length.
  //
  // Each carries the name its volume takes, since the layer holds both
  // poly and copper and they are not numbered into one series.
  struct BottomWallPrism {
    std::string         name;
    std::string         material;
    std::vector<double> uVerts;
    std::vector<double> vVerts;
    double              length;        // the sweep, its thickness in y
    std::string         orientation;
    // Placed by an ANCHOR on the cap's origin, not by a center. The
    // sweep is centered on that point, as G4ExtrudedSolid centers an
    // extrusion on its placement point, so the anchor's y is the
    // mid-plane of the thickness rather than a face.
    CLHEP::Hep3Vector   anchor;
  };

  class BottomWall {
  public:

    BottomWall(bool build,
               BottomWallPlate const & basePlate,
               std::vector<BrickWall> const & leadLayers,
               std::vector<BottomWallPrism> const & prisms
               ) :
      _build(build),
      _basePlate(basePlate),
      _leadLayers(leadLayers),
      _prisms(prisms)
    {
    }

    bool build() const {return _build;}

    // The steel plate the house stands on.
    BottomWallPlate const & basePlate() const {return _basePlate;}

    // The two lead layers, in the order they are stacked. Each holds
    // its course bricks and its strays together.
    std::vector<BrickWall> const & leadLayers() const {return _leadLayers;}

    // The swept pieces, in the order the layers are built: the lower
    // poly L, the upper poly L, then the two copper prisms.
    std::vector<BottomWallPrism> const & prisms() const {return _prisms;}

    // Genreflex can't do persistency of vector<BottomWall> without a
    // default constructor
    BottomWall() {}

  private:

    bool _build;

    BottomWallPlate              _basePlate;
    std::vector<BrickWall>       _leadLayers;
    std::vector<BottomWallPrism> _prisms;
  };

}

#endif/*STMGeom_BottomWall_hh*/
