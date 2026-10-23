#ifndef STMGeom_FrontShieldingRight_hh
#define STMGeom_FrontShieldingRight_hh

// The right half of the STM front shielding: the section behind the
// spot-size collimator.
//
// Five layers stacked into the wall -- lead, borated poly, lead,
// borated poly, lead -- standing on four aluminium blocks, with two
// collimator pipes and a copper plate in front.
//
// The lead layers are BrickWalls, each a corner and three directions
// with its bricks already placed. The poly sheets are held here
// rather than in BrickWall because the interleaving is this section's
// own: another wall splits its sheets differently, or not at all, so
// there is no shared sequence to factor out.
//
// Every position arrives resolved from STMMaker. The section is fixed
// by brickEndX, the bore height above the baseplate, and the gap
// behind the cradle; everything else follows from butting pieces
// together, so no placement here was measured off the CAD.
//
// Author: Yongyi Wu

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/STMGeom/inc/BrickWall.hh"

namespace mu2e {

  // A poly plate between the lead layers. Each is listed on its own,
  // with its own size and bores.
  struct FrontShieldingRightSheet {
    std::string       material;
    CLHEP::Hep3Vector halfDim;
    CLHEP::Hep3Vector center;
    std::string       orientation;
    std::vector<int>  bores;           // bore ids, indexing the layer's bores
    std::vector<double> boreRadius;    // one per entry in bores
  };

  // The copper plate in front: a trapezoid swept along one axis and
  // bored on the beam axes. Its outline is given in the
  // ExtShieldDownstream form, a u/v polygon plus a sweep length.
  struct FrontShieldingRightPlate {
    std::string         material;
    std::vector<double> uVerts;
    std::vector<double> vVerts;
    double              length;
    std::string         orientation;
    // Placed by an ANCHOR, not a center: the outline is drawn in the
    // third quadrant so its local origin is the +u end of the long
    // side, and that corner sits a stated distance from the brick end,
    // whereas a trapezoid's center lines up with nothing.
    CLHEP::Hep3Vector   anchor;
    // Which bores pass through the plate. The ids matter rather than
    // the order, since a bore may be moved off the beam axis by its own
    // offset and the construction code must read that center.
    std::vector<int>    bores;
    std::vector<double> holeRadius;    // one per entry in bores
  };

  // A collimator pipe on one of the beam axes.
  struct FrontShieldingRightPipe {
    std::string       material;
    double            rIn;
    double            rOut;
    double            halfLength;
    CLHEP::Hep3Vector center;
  };

  class FrontShieldingRight {
  public:

    FrontShieldingRight(bool build,
                        std::vector<BrickWall> const & leadLayers,
                        std::vector<FrontShieldingRightSheet> const & sheets,
                        std::string const & blockMaterial,
                        CLHEP::Hep3Vector const & blockHalfDim,
                        std::vector<CLHEP::Hep3Vector> const & blockCenter,
                        std::vector<FrontShieldingRightPipe> const & pipes,
                        FrontShieldingRightPlate const & plate,
                        double backZ
                        ) :
      _build(build),
      _leadLayers(leadLayers),
      _sheets(sheets),
      _blockMaterial(blockMaterial),
      _blockHalfDim(blockHalfDim),
      _blockCenter(blockCenter),
      _pipes(pipes),
      _plate(plate),
      _backZ(backZ)
    {
    }

    bool build() const {return _build;}

    // The three lead layers, in depth order.
    std::vector<BrickWall> const & leadLayers() const {return _leadLayers;}

    // The poly plates, in depth order.
    std::vector<FrontShieldingRightSheet> const & sheets() const {return _sheets;}

    // The aluminium blocks the wall stands on. They span the same
    // width as the courses and their undersides land on the
    // baseplate, so their centers are derived, not measured.
    std::string const & blockMaterial() const {return _blockMaterial;}
    CLHEP::Hep3Vector const & blockHalfDim() const {return _blockHalfDim;}
    std::vector<CLHEP::Hep3Vector> const & blockCenter() const {return _blockCenter;}

    std::vector<FrontShieldingRightPipe> const & pipes() const {return _pipes;}
    FrontShieldingRightPlate const & plate() const {return _plate;}

    // The z of the section's downstream face -- the back of the copper
    // lining. Sections further down the house position themselves
    // against this, so it is resolved once in STMMaker, where the
    // lining's rotation is in hand, rather than re-derived from an
    // outline whose extents depend on how it is turned.
    double backZ() const {return _backZ;}

    // Genreflex can't do persistency of vector<FrontShieldingRight>
    // without a default constructor
    FrontShieldingRight() {}

  private:

    bool _build;

    std::vector<BrickWall> _leadLayers;
    std::vector<FrontShieldingRightSheet> _sheets;

    std::string       _blockMaterial;
    CLHEP::Hep3Vector _blockHalfDim;
    std::vector<CLHEP::Hep3Vector> _blockCenter;

    std::vector<FrontShieldingRightPipe> _pipes;
    FrontShieldingRightPlate             _plate;
    double                               _backZ;
  };

}

#endif/*STMGeom_FrontShieldingRight_hh*/
