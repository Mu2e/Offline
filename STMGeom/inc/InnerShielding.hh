#ifndef STMGeom_InnerShielding_hh
#define STMGeom_InnerShielding_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023
// Update: Yongyi Wu Sept-Oct 2026, updated (hand-stacked) inner shielding

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/STMGeom/inc/BrickWall.hh"

namespace mu2e {

  // ---- the updated, hand-stacked inner shielding ---------------------
  //
  // What lines the cavity the detectors sit in: six copper pieces and
  // eleven lead ones. Unlike the walls this is not a stack of layers --
  // every piece is placed on its own -- so there is no layer sequence,
  // just the three lists below.
  //
  // Some pieces are bored on the LaBr beam axis, naming the bores they
  // take by number as everything else in the house does. The bore
  // follows the collimator, so STMMaker resolves its center from
  // offset_Spot.

  // A prism swept along one axis, given in the ExtShieldDownstream
  // form: a u/v cap polygon plus a sweep length.
  //
  // Placed by an ANCHOR on its cap origin, since the outlines are
  // traced from the vertex nearest the section reference. The sweep is
  // centered on that point, as G4ExtrudedSolid centers an extrusion on
  // its placement point.
  struct InnerShieldingPrism {
    std::string         name;
    std::string         material;
    std::vector<double> uVerts;
    std::vector<double> vVerts;
    double              length;        // the sweep
    std::string         orientation;
    CLHEP::Hep3Vector   anchor;
    std::vector<int>    bores;         // bore ids, indexing bores()
  };

  // A plain box: the two copper plates, and the three lead pieces that
  // are turned 45 degrees about y and so cannot be stated as leadBrick
  // types.
  struct InnerShieldingBox {
    std::string       name;
    std::string       material;
    CLHEP::Hep3Vector halfDim;
    CLHEP::Hep3Vector center;
    std::string       orientation;
    std::vector<int>  bores;           // bore ids, indexing bores()
    // Set where a nominal size would touch its neighbours once turned.
    // The construction code backs each half-dimension off, so halfDim
    // above stays the nominal value the config states.
    bool              nudge;
  };

  // The standard lead bricks are BrickWallBrick, as in every wall: this
  // section places each one individually rather than expanding courses,
  // but a placed brick is a placed brick, and sharing the type lets the
  // construction code use the same helper.

  class InnerShielding {
  public:

    // The earlier description. Its geometry is written out in
    // constructSTM.cc, so the class carries only the switch.
    InnerShielding(bool build):
      _build(build)
    {}

    // The updated, hand-stacked description. A separate constructor
    // rather than defaulted arguments, so neither can be built with a
    // quantity that does not apply to it.
    //
    // Every position arrives resolved from STMMaker, fixed by one
    // reference -- the bottom wall's, moved in 0.5 in and up 6 in.
    InnerShielding(bool build,
                   std::vector<InnerShieldingPrism> const & prisms,
                   std::vector<InnerShieldingBox> const & boxes,
                   std::vector<BrickWallBrick> const & bricks,
                   std::vector<BrickWallBore> const & bores):
      _build(build),
      _prisms(prisms),
      _boxes(boxes),
      _bricks(bricks),
      _bores(bores)
    {}

    bool   build()                               const { return _build;  }

    // ---- updated (hand-stacked) only; empty for the earlier one -----

    // The swept pieces: four copper, and on the lead side the angled
    // triangle, the bored trapezoid and the two clipped prisms.
    std::vector<InnerShieldingPrism> const & prisms() const {return _prisms;}

    // The boxes: two copper plates and the three angled lead pieces.
    std::vector<InnerShieldingBox> const & boxes() const {return _boxes;}

    // The standard lead bricks placed individually.
    std::vector<BrickWallBrick> const & bricks() const {return _bricks;}

    // The holes through this section, named by number from 1.
    std::vector<BrickWallBore> const & bores() const {return _bores;}

    InnerShielding() {}
  private:

    bool               _build;

    std::vector<InnerShieldingPrism> _prisms;
    std::vector<InnerShieldingBox>   _boxes;
    std::vector<BrickWallBrick>      _bricks;
    std::vector<BrickWallBore>       _bores;

  };

}

#endif/*STMGeom_InnerShielding_hh*/
