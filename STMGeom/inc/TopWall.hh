#ifndef STMGeom_TopWall_hh
#define STMGeom_TopWall_hh

// The roof of the STM shield house.
//
// Three layers worked downward in -y from the top of the stack: a pair
// of borated poly sheets, a layer of lead bricks laid flat, and a pair
// of aluminium plates. 29 pieces, none of them bored.
//
// The lead layer is a BrickWall like the side walls', but laid flat:
// its courses run side by side in x and each course runs along z, so
// the 2 in of a brick is the layer's thickness in y. Two of the five
// courses start 3 in along z from the rest, which is why each carries
// its own offset rather than sharing one origin.
//
// Every position arrives resolved from STMMaker. The wall is fixed by
// one reference -- inboard of the right wall in x, at that wall's
// reference z, and at the top of the stack in y -- so each piece
// follows from butting against its neighbour.
//
// Author: Yongyi Wu

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/STMGeom/inc/BrickWall.hh"

namespace mu2e {

  // Its flat plates are WallSheet, shared with the other walls.

  class TopWall {
  public:

    TopWall(bool build,
            std::vector<BrickWall> const & leadLayers,
            std::vector<WallSheet> const & sheets
            ) :
      _build(build),
      _leadLayers(leadLayers),
      _sheets(sheets)
    {
    }

    bool build() const {return _build;}

    // The lead. One layer, held as a vector so the construction code
    // walks it the same way it walks the side walls'.
    std::vector<BrickWall> const & leadLayers() const {return _leadLayers;}

    // The flat pieces, in the order the layers are built: the two poly
    // sheets, then the two aluminium plates.
    std::vector<WallSheet> const & sheets() const {return _sheets;}

    // Genreflex can't do persistency of vector<TopWall> without a
    // default constructor
    TopWall() {}

  private:

    bool _build;

    std::vector<BrickWall> _leadLayers;
    std::vector<WallSheet> _sheets;
  };

}

#endif/*STMGeom_TopWall_hh*/
