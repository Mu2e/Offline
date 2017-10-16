#ifndef Alignment_AlignmentObj_HH
#define Alignment_AlignmentObj_HH
// Alignment/inc/AlignmentObj.hh
// This is the header file for the AlignmentObj class, which 
// represents the characterization of the misalignment
// of a geometry object.  Objects of this class not only 
// represent the misalignments observed in the as-built geometry,
// but can also be used either to create simulated data with
// the same misalignments or to correct data for the misalignments.

// Original Author:  David Norvil Brown (Louisville), September 2017

// In the first version, we assume all misalignments can be characterized 
// as a combination of a rotation (relative to its own axes, not external
// axes) and a displacement.  Note that this works for placement error
// misalignment but not deformation misalignment.

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "Alignment/inc/ShapeDetail.hh"

namespace mu2e {

  class AlignmentObj {
  public:
    AlignmentObj() {_displacement.set(0,0,0); _rotation=CLHEP::HepRotation::IDENTITY; _detail=0;}
    AlignmentObj( CLHEP::Hep3Vector & disp, CLHEP::HepRotation & rot, 
		  ShapeDetail* deet = 0 ) :
      _displacement(disp),
      _rotation(rot),
      _detail( deet ){}

    CLHEP::Hep3Vector  displacement() const { return _displacement;}
    CLHEP::HepRotation rotation()     const { return _rotation;    }
    ShapeDetail*       detail()       const { return _detail;      }

    void  setDetail( ShapeDetail* theDeets ) { _detail = theDeets; }
    void  setDisplacement ( CLHEP::Hep3Vector & aDisp ) { _displacement = aDisp;}
    void  setRotation ( CLHEP::HepRotation & aRot ) { _rotation = aRot; }

  private:
    CLHEP::Hep3Vector  _displacement;
    CLHEP::HepRotation _rotation;
    ShapeDetail*       _detail;
  };
} // end of namespace mu2e
#endif  //  Alignment_AlignmentObj_HH
