#ifndef STMGeom_LeadBrick_hh
#define STMGeom_LeadBrick_hh

// The standard lead bricks the STM shield house is stacked from.
//
// A shared primitive rather than a property of any one structure: the
// same sizes recur throughout the shielding, so their dimensions,
// material and wear live here.
//
// The sizes are a numbered list rather than named members, so a course
// anywhere in the house is a list of type numbers and a new size is a
// new entry rather than a change to every caller. A type is a SHAPE
// only -- no orientation, since the same brick lies flat in one wall
// and stands on end in another, and no bores, since which holes pass
// through a brick depends on where it sits.
//
// Wear is taken off each face, so a brick shrinks while staying
// centered where it was placed: a stack keeps its nominal pitch and
// the wear opens as gaps. The as-delivered sizes are exact imperial.
//
// The wear is not isotropic -- one value for whichever of the brick's
// own axes ends up vertical, another for the other two -- so worn()
// takes the placement rotation and the construction code needs a solid
// per (type, orientation) pair.
//
// Author: Yongyi Wu

#include <cmath>
#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class LeadBrick {
  public:

    LeadBrick(std::vector<std::string> const & names,
              std::vector<CLHEP::Hep3Vector> const & dims,
              double wearY,
              double wearXZ,
              std::string const & material
              ) :
      _names(names),
      _dims(dims),
      _wearY(wearY),
      _wearXZ(wearXZ),
      _material(material)
    {
    }

    // How many sizes are defined. Types are numbered from 1.
    size_t nTypes() const {return _dims.size();}

    // The name a type was given in the geometry file, e.g. "2x4x16".
    std::string const & name(int type) const {return _names.at(type-1);}

    // Nominal outside dimensions of a type, before wear.
    CLHEP::Hep3Vector const & dim(int type) const {return _dims.at(type-1);}

    // The same with the wear taken off each face, in the brick's OWN
    // frame. This is what the G4Box should be built from.
    //
    // rot is the rotation handed to G4PVPlacement, which takes it as the
    // frame rotation, so the brick itself is turned by rot^-1 and its
    // local axis k lands on rot^-1 * e_k in Mu2e. The local axis that
    // lands closest to Mu2e y takes wearY and the other two take wearXZ.
    // For the 90 degree turns used here each local axis lands exactly on
    // a Mu2e axis, so that is unambiguous.
    CLHEP::Hep3Vector worn(int type, CLHEP::HepRotation const & rot) const {
      CLHEP::Hep3Vector const & d = _dims.at(type-1);
      double w[3];
      for (int k = 0; k < 3; ++k) {
        const CLHEP::Hep3Vector axis =
          rot.inverse() * CLHEP::Hep3Vector(k == 0, k == 1, k == 2);
        w[k] = (std::abs(axis.y()) > 0.5) ? _wearY : _wearXZ;
      }
      return CLHEP::Hep3Vector(d.x() - 2.*w[0],
                               d.y() - 2.*w[1],
                               d.z() - 2.*w[2]);
    }

    double              wearY()    const {return _wearY;}
    double              wearXZ()   const {return _wearXZ;}
    std::string const & material() const {return _material;}

    // Genreflex can't do persistency of vector<LeadBrick> without a
    // default constructor
    LeadBrick() {}

  private:

    std::vector<std::string>       _names;
    std::vector<CLHEP::Hep3Vector> _dims;
    double                         _wearY;
    double                         _wearXZ;
    std::string                    _material;
  };

}

#endif/*STMGeom_LeadBrick_hh*/
