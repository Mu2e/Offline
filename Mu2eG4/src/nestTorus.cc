//
// Free function to create and place a new G4Torus, place inside a logical volume.
// 
// $Id: nestTorus.cc,v 1.1.1.1
// $Author: A Chandra
// $Date: 2010/03/15
//

#include <string>

#include "Mu2eG4/inc/nestTorus.hh"

#include "G4Torus.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

using namespace std;

namespace mu2e {
 
  //
  // Create and place a G4Torus inside a logical volume.
  // 
  VolumeInfo nestTorus ( string const& name,
                        double param[5],
                        G4Material* material,
                        G4RotationMatrix* rot,
                        G4ThreeVector const& offset,
                        G4LogicalVolume* parent,
                        int copyNo,
                        G4Colour color,
                        bool forceSolid
                        ){
    
    VolumeInfo info;
    
    info.solid   = new G4Torus( name, param[0], param[1], param[2], param[3], param[4]  );
    
    info.logical = new G4LogicalVolume( info.solid, material, name); 
    
    info.physical =  new G4PVPlacement( rot, offset, info.logical, name, parent, 0, copyNo);
    
    G4VisAttributes* visAtt = new G4VisAttributes(true, color);
    visAtt->SetForceSolid(forceSolid);

    visAtt->SetForceAuxEdgeVisible (false);

    info.logical->SetVisAttributes(visAtt);

    return info;
  }

}
