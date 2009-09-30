//
// Free function to create a new G4 Box, placed inside a logical volume.
// 
// $Id: nestBox.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <string>

#include "Mu2eG4/inc/nestBox.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

using namespace std;

namespace mu2e {
 
  //
  // Place a box inside a logical volume.
  // 
  VolumeInfo nestBox ( string const& name,
		       double halfDim[3],
		       G4Material* material,
		       G4RotationMatrix* rot,
		       G4ThreeVector const& offset,
		       G4LogicalVolume* parent,
		       int copyNo,
		       G4Colour color,
		       bool forceSolid
		       ){
    
    VolumeInfo info;
    
    info.solid   = new G4Box( name, halfDim[0], halfDim[1], halfDim[2] );
    
    info.logical = new G4LogicalVolume( info.solid, material, name); 
    
    info.physical =  new G4PVPlacement( rot, offset, info.logical, name, parent, 0, copyNo);
    
    G4VisAttributes* visAtt = new G4VisAttributes(true, color);
    visAtt->SetForceSolid(forceSolid);
    info.logical->SetVisAttributes(visAtt);
    
    return info;
  }

}
