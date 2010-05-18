//
// Free function to create and place a new G4Tubs, place inside a logical volume.
// 
// $Id: nestTubs.cc,v 1.3 2010/05/18 21:16:28 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/05/18 21:16:28 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) See note in Mu2eWorld about how G4VisAttributes leaks memory.
//    Should address this sometime.

#include <string>

#include "Mu2eG4/inc/nestTubs.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

using namespace std;

namespace mu2e {
 
  //
  // Create and place a G4Tubs inside a logical volume.
  // 
  VolumeInfo nestTubs ( string const& name,
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
    
    info.solid   = new G4Tubs( name, param[0], param[1], param[2], param[3], param[4]  );
    
    info.logical = new G4LogicalVolume( info.solid, material, name); 
    
    info.physical =  new G4PVPlacement( rot, offset, info.logical, name, parent, 0, copyNo);

    // This leaks. See note 1.
    G4VisAttributes* visAtt = new G4VisAttributes(true, color);
    visAtt->SetForceSolid(forceSolid);

    // If I do not do this, then the rendering depends on what happens in
    // other parts of the code;  is there a G4 bug that causes something to be
    // unitialized?
    visAtt->SetForceAuxEdgeVisible (false);

    info.logical->SetVisAttributes(visAtt);

    return info;
  }

}
