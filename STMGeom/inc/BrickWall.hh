#ifndef STMGeom_BrickWall_hh
#define STMGeom_BrickWall_hh

// A layer in the wall of the STM shield house, stacked from standard
// lead bricks, and can interleave with other sheet layers.
//
// Shared description for most of the house. A wall layer is fixed
// by a corner and three directions:
//        ____________________________
// pitch /     /       /             /|
//   ^  /_____/_______/_____________/ |
//   |  |_____|_______|_____________|/|
//   |  |____________________|______|/
//  origin       course-->
//
//   origin      the corner every course is flush against
//   courseDir   the way a course runs from that corner
//   pitchDir    the way courses stack across the layer
//   depthDir    which side of the origin the layer occupies
//
// Carrying the directions as vectors rather than axis names settles
// the signs: the front wall's first course is origin (brickEndX, bore,
// front face) with courseDir -x, pitchDir +y, depthDir +z; a right
// wall is the same class with courseDir +z and depthDir -x.
//
// This is ONE layer of bricks. A wall is several of these interleaved
// with sheets, and the section owns that sequence, since sheets differ
// from layer to layer.
//
// A brick's place follows from butting its course's type list end to
// end against the origin, so no position here is a measured one.
// STMMaker does that expansion and hands over finished placements.
//
// Bores are kept apart from brick types on purpose: a type is a shape
// only, so a brick catching one hole and a brick catching two remain
// the same type. The bore map is given rather than computed so that it
// can be checked against the geometry -- if the origin moves far
// enough to slide a hole onto a neighbouring brick, the mismatch is
// caught instead of quietly boring the wrong piece.
//
// Author: Yongyi Wu

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  // One brick, as placed.
  struct BrickWallBrick {
    int                 type;          // index into LeadBrick's types
    CLHEP::Hep3Vector   center;
    std::string         orientation;   // OrientationResolver code
    std::vector<int>    bores;         // bore ids through this brick, may be empty
  };

  // One flat plate in a wall: the poly and copper sheets and the
  // aluminium plates. Shared by all four walls, which differ in which
  // plates they have and where those sit, never in what a plate is.
  //
  // No bores and no orientation -- halfDim is already in the Mu2e
  // frame, so which of the three is the thickness is read off it.
  // STMMaker sets each name, where which-plate-is-which is known;
  // numbering them in the construction code would put poly, copper and
  // aluminium in one series and lose the material prefix.
  struct WallSheet {
    std::string       name;
    std::string       material;
    CLHEP::Hep3Vector halfDim;
    CLHEP::Hep3Vector center;
  };

  // One hole. It is pinned to a beam axis rather than to the wall, so
  // it keeps its place when the wall around it is retuned.
  struct BrickWallBore {
    std::string       axis;            // which beam axis it follows
    double            radius;
    CLHEP::Hep3Vector offset;          // correction from coaxial
    CLHEP::Hep3Vector center;          // resolved by STMMaker
  };

  class BrickWall {
  public:

    BrickWall(bool build,
              CLHEP::Hep3Vector const & origin,
              CLHEP::Hep3Vector const & courseDir,
              CLHEP::Hep3Vector const & pitchDir,
              CLHEP::Hep3Vector const & depthDir,
              std::vector<BrickWallBrick> const & bricks,
              std::vector<BrickWallBore>  const & bores
              ) :
      _build(build),
      _origin(origin),
      _courseDir(courseDir),
      _pitchDir(pitchDir),
      _depthDir(depthDir),
      _bricks(bricks),
      _bores(bores)
    {
    }

    bool build() const {return _build;}

    // The frame the wall was built in. Kept so that the construction
    // code and any later check can say where a piece sits without
    // guessing which plane the wall lies in.
    CLHEP::Hep3Vector const & origin()    const {return _origin;}
    CLHEP::Hep3Vector const & courseDir() const {return _courseDir;}
    CLHEP::Hep3Vector const & pitchDir()  const {return _pitchDir;}
    CLHEP::Hep3Vector const & depthDir()  const {return _depthDir;}

    std::vector<BrickWallBrick> const & bricks() const {return _bricks;}
    std::vector<BrickWallBore>  const & bores()  const {return _bores;}

    // Genreflex can't do persistency of vector<BrickWall> without a
    // default constructor
    BrickWall() {}

  private:

    bool _build;

    CLHEP::Hep3Vector _origin;
    CLHEP::Hep3Vector _courseDir;
    CLHEP::Hep3Vector _pitchDir;
    CLHEP::Hep3Vector _depthDir;

    std::vector<BrickWallBrick> _bricks;
    std::vector<BrickWallBore>  _bores;
  };

}

#endif/*STMGeom_BrickWall_hh*/
