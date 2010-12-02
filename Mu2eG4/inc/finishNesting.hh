#ifndef FINISHNESTING_HH
#define FINISHNESTING_HH
//
// Free function to be used by the nest... functions
// 
// $Id: finishNesting.hh,v 1.1 2010/12/02 17:43:54 genser Exp $
// $Author: genser $ 
// $Date: 2010/12/02 17:43:54 $
//
// Original author KLG
//


//class G4RotationMatrix;
//class G4ThreeVector;
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class VolumeInfo;

class G4Material;
class G4LogicalVolume;
class G4Colour;

namespace mu2e {

  void finishNesting(VolumeInfo& info,
                     G4Material* material,
                     G4RotationMatrix* rot,
                     G4ThreeVector const & offset,
                     G4LogicalVolume* parent,
                     int copyNo,
                     bool const isVisible,
                     G4Colour const color,
                     bool const forceSolid,
                     bool const forceAuxEdgeVisible,
                     bool const placePV,
                     bool const doSurfCheck
                     );

}

#endif
