#ifndef NESTTORUS_HH
#define NESTTORUS_HH
//
// Free function to create and place a new G4Torus, place inside a logical volume.
// 
// $Id: v 1.1.1.1
// $Author: genser $ 
// $Date: 2010/03/15
//

#include <string>
#include <vector>

#include "Mu2eG4/inc/VolumeInfo.hh"

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4CSGSolid;

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"


namespace mu2e {

  VolumeInfo nestTorus ( std::string const& name,
                         double halfDim[5],
                         G4Material* material,
                         G4RotationMatrix* rot,
                         G4ThreeVector const& offset,
                         G4LogicalVolume* parent,
                         int copyNo,
                         G4Colour color,
                         bool forceSolid,
                         bool doSurfaceCheck
                         );
  


  // Alternate argument list, using a vector for the half dimensions.
  //
  inline VolumeInfo nestTorus ( std::string const& name,
                                std::vector<double>&  halfDim,
                                G4Material* material,
                                G4RotationMatrix* rot,
                                G4ThreeVector& offset,
                                G4LogicalVolume* parent,
                                int copyNo,
                                G4Colour color,
                                bool forceSolid,
                                bool doSurfaceCheck
                                ){
    return nestTorus( name, 
                      &halfDim[0],
                      material,
                      rot,
                      offset,
                      parent,
                      copyNo,
                      color,
                      forceSolid,
                      doSurfaceCheck
                      );
  }
}

#endif
