#ifndef STMGeom_RightWall_hh
#define STMGeom_RightWall_hh

// The right wall of the STM shield house: the wall on the +x side,
// running downstream from the front shielding.
//
// Five layers worked inward along -x from one reference vertex -- an
// outer poly pair, a lead layer, an inner poly pair, a second lead
// layer, and a copper sheet -- 26 pieces. Nothing here is bored.
//
// The lead layers are BrickWalls, as in the front shielding, with
// courses running along z and the layers stacking along -x.
//
// Every position arrives resolved from STMMaker. The wall is fixed by
// the reference vertex -- flush with the +z side of the front
// shielding's last lead layer, at that layer's +x, -y corner -- and by
// the baseplate, so each piece follows from butting against its
// neighbour rather than from a measurement of its own.
//
// Author: Yongyi Wu

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/STMGeom/inc/BrickWall.hh"

namespace mu2e {

  // Its flat plates are WallSheet, shared with the other walls.

  // The L-shaped edge piece at the upstream end: an outline in the
  // y-z plane swept along x, given in the ExtShieldDownstream form.
  struct RightWallEdgePrism {
    std::string         material;
    std::vector<double> uVerts;
    std::vector<double> vVerts;
    double              length;        // the sweep, its thickness in x
    std::string         orientation;
    // Placed by an ANCHOR, not a center: the outline is drawn from the
    // corner where it meets the front shielding, and an L has no center
    // that lines up with anything. The sweep is centered on the anchor,
    // as G4ExtrudedSolid centers an extrusion on its placement point.
    CLHEP::Hep3Vector   anchor;
  };

  class RightWall {
  public:

    RightWall(bool build,
              std::vector<BrickWall> const & leadLayers,
              std::vector<WallSheet> const & sheets,
              RightWallEdgePrism const & edgePrism
              ) :
      _build(build),
      _leadLayers(leadLayers),
      _sheets(sheets),
      _edgePrism(edgePrism)
    {
    }

    bool build() const {return _build;}

    // The two lead layers, worked inward: four courses of three, then
    // three courses of three.
    std::vector<BrickWall> const & leadLayers() const {return _leadLayers;}

    // The flat sheets, in the order the layers are built: the outer top
    // sheet, the inner longwall and top sheet, then the copper sheet.
    std::vector<WallSheet> const & sheets() const {return _sheets;}

    RightWallEdgePrism const & edgePrism() const {return _edgePrism;}

    // Genreflex can't do persistency of vector<RightWall> without a
    // default constructor
    RightWall() {}

  private:

    bool _build;

    std::vector<BrickWall> _leadLayers;
    std::vector<WallSheet> _sheets;
    RightWallEdgePrism     _edgePrism;
  };

}

#endif/*STMGeom_RightWall_hh*/
