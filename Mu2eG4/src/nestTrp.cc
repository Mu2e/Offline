//
// Free function to create a new G4 Trp, placed inside a logical volume.
// 
// $Id: nestTrp.cc,v 1.2 2010/07/30 19:43:25 genser Exp $
// $Author: genser $ 
// $Date: 2010/07/30 19:43:25 $
//
// Original author Krzysztof Genser based on Rob Kutschke' nestBox
//

#include <string>

#include "Mu2eG4/inc/nestTrp.hh"

//#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

using namespace std;

namespace mu2e {
 
  //
  // Place a trp inside a logical volume.
  // 
  VolumeInfo nestTrp ( string const& name,
                       double const halfDim[5],
                       G4Material* material,
                       G4RotationMatrix* rot,
                       G4ThreeVector const& offset,
                       G4LogicalVolume* parent,
                       int copyNo,
                       G4Colour color,
                       bool forceSolid
                       ){
    
    VolumeInfo info;

//    info.solid   = new G4Box( name, halfDim[0], halfDim[1], halfDim[2] );
//    z,y,x, smallerx

    info.solid   = new G4Trd ( name,
			       halfDim[4], 
			       halfDim[3], 
			       halfDim[2],		       
			       halfDim[2],
			       halfDim[1]
			       );

    info.logical  = new G4LogicalVolume( info.solid, material, name); 
    
    info.physical = new G4PVPlacement( rot, offset, info.logical, name, parent, 0, copyNo);
    
    G4VisAttributes* visAtt = new G4VisAttributes(true, color);
    visAtt->SetForceSolid(forceSolid);
    visAtt->SetForceAuxEdgeVisible(false);
    info.logical->SetVisAttributes(visAtt);
    
    return info;
  }

}
