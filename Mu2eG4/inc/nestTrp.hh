#ifndef Mu2eG4_nestTrp_hh
#define Mu2eG4_nestTrp_hh
//
// Free function to create a new G4 Trp, placed inside a logical volume.
//
//
// Original author Krzysztof Genser based on Rob Kutschke's nestBox
//

#include <string>
#include <vector>

#include "Mu2eG4Helper/inc/VolumeInfo.hh"

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4CSGSolid;

#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4Colour.hh"

namespace mu2e {

  VolumeInfo nestTrp ( std::string const& name,
                       double const halfDim[5],
                       G4Material* material,
                       G4RotationMatrix const* rot,
                       G4ThreeVector const& offset,
                       G4LogicalVolume* parent,
                       int copyNo,
                       bool const isVisible,
                       G4Colour const color,
                       bool const forceSolid,
                       bool const forceAuxEdgeVisible,
                       bool const placePV,
                       bool const doSurfaceCheck
                       );



  // Alternate argument list, using a vector for the half dimensions.
  inline VolumeInfo nestTrp ( std::string const& name,
                              std::vector<double> const&  halfDim,
                              G4Material* material,
                              G4RotationMatrix const* rot,
                              G4ThreeVector const& offset,
                              G4LogicalVolume* parent,
                              int copyNo,
                              bool const isVisible,
                              G4Colour const color,
                              bool const forceSolid,
                              bool const forceAuxEdgeVisible,
                              bool const placePV,
                              bool const doSurfaceCheck
                              ){
    return nestTrp( name,
                    &halfDim[0],
                    material,
                    rot,
                    offset,
                    parent,
                    copyNo,
                    isVisible,
                    color,
                    forceSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck
                    );
  }
}

#endif /* Mu2eG4_nestTrp_hh */
