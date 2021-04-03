//
// Free function to create a new G4 Trp, placed inside a logical volume.
//
//
// Original author Krzysztof Genser based on Rob Kutschke's nestBox
//

#include <string>

#include "Mu2eG4/inc/nestTrp.hh"
#include "Mu2eG4/inc/finishNesting.hh"

#include "Geant4/G4Trd.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4ThreeVector.hh"

using namespace std;

namespace mu2e {

  //
  // Place a trp inside a logical volume.
  //
  VolumeInfo nestTrp ( string const& name,
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
                       ){

    VolumeInfo info;

    info.name     = name;

//    info.solid   = new G4Box( name, halfDim[0], halfDim[1], halfDim[2] );
//    z,y,x, smallerx

     info.solid   = new G4Trd ( name,
                                halfDim[4],
                                halfDim[3],
                                halfDim[2],
                                halfDim[2],
                                halfDim[1]
                                );

     finishNesting(info,
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

     return info;
  }

}
