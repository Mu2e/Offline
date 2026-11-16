#ifndef STMGeom_LeftWall_hh
#define STMGeom_LeftWall_hh

// The left wall of the STM shield house: the wall along -x, running
// downstream from the front shielding.
//
// The mirror of the right wall, with the same five kinds of layer and
// the same sheet shapes, but no L-shaped edge piece: 26 pieces. Nothing
// here is bored.
//
// Layers are listed outside to inside and walk +x from the reference,
// which is the -x face of the outermost layer, so every piece is on the
// same side of it. The lead layers are BrickWalls with courses running
// along z, as in the right wall.
//
// Every position arrives resolved from STMMaker. The wall is fixed by
// that reference and by the baseplate, so each piece follows from
// butting against its neighbour rather than from a measurement of its
// own.
//
// Deliberately its own class rather than a generalisation of RightWall:
// the four walls hold different sets of parts -- this one has no edge
// prism, the floor has a baseplate and prisms -- so one class would
// carry empty members for most of them. What they do share is factored
// out instead: the lead layers through BrickWall and layLayer, the
// plates through WallSheet.
//
// Author: Yongyi Wu

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/STMGeom/inc/BrickWall.hh"

namespace mu2e {

  // Its flat plates are WallSheet, shared with the other walls.

  class LeftWall {
  public:

    LeftWall(bool build,
             std::vector<BrickWall> const & leadLayers,
             std::vector<WallSheet> const & sheets
             ) :
      _build(build),
      _leadLayers(leadLayers),
      _sheets(sheets)
    {
    }

    bool build() const {return _build;}

    // The two lead layers, outside to inside: four courses of three,
    // then three courses of three.
    std::vector<BrickWall> const & leadLayers() const {return _leadLayers;}

    // The flat sheets, in the order the layers are built: the outer
    // longwall and top edge, the inner pair, then the copper sheet.
    std::vector<WallSheet> const & sheets() const {return _sheets;}

    // Genreflex can't do persistency of vector<LeftWall> without a
    // default constructor
    LeftWall() {}

  private:

    bool _build;

    std::vector<BrickWall> _leadLayers;
    std::vector<WallSheet> _sheets;
  };

}

#endif/*STMGeom_LeftWall_hh*/
