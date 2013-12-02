#ifndef Mu2eG4_finishNesting_hh
#define Mu2eG4_finishNesting_hh
//
// Free function to be used by the nest... functions
//
// $Id: finishNesting.hh,v 1.6 2013/12/02 20:06:13 genser Exp $
// $Author: genser $
// $Date: 2013/12/02 20:06:13 $
//
// Original author KLG
//

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4Material;
class G4LogicalVolume;
class G4Colour;

namespace mu2e {

  class VolumeInfo;

  void finishNesting(VolumeInfo& info,
                     G4Material* material,
                     G4RotationMatrix const* rot,
                     G4ThreeVector const & offset,
                     G4LogicalVolume* parent,
                     int copyNo,
                     bool const isVisible,
                     G4Colour const color,
                     bool const forceSolid,
                     bool const forceAuxEdgeVisible,
                     bool const placePV,
                     bool const doSurfaceCheck,
                     bool const verbose = false
                     );

}

#endif /* Mu2eG4_finishNesting_hh */
