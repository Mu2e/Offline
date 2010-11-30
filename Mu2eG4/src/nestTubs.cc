//
// Free function to create and place a new G4Tubs, place inside a logical volume.
// 
// $Id: nestTubs.cc,v 1.6 2010/11/30 16:38:24 genser Exp $
// $Author: genser $ 
// $Date: 2010/11/30 16:38:24 $
//
// Original author Rob Kutschke
//

#include <string>

// Mu2e includes
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/G4Helper.hh"

// G4 includes
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
                        bool isVisible,
                        G4Colour color,
                        bool forceSolid,
                        bool forceAuxEdgeVisible,
                        bool placePV,
                        bool doSurfCheck
                        ){
    

    G4Helper    & _helper = *(edm::Service<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    VolumeInfo info;
    
    info.solid    = new G4Tubs( name, param[0], param[1], param[2], param[3], param[4]  );
    
    info.logical  = new G4LogicalVolume( info.solid, material, name); 

    info.physical  =  placePV ?
      info.physical = new G4PVPlacement( rot, offset, info.logical, name, parent, 0, copyNo, 
                                         doSurfCheck)
      :
      0;

    info.name     = name;

    if (!isVisible) {

      info.logical->SetVisAttributes(G4VisAttributes::Invisible);

    } else {
      
      G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, color));
      visAtt->SetForceSolid(forceSolid);
      // If I do not do this, then the rendering depends on what happens in
      // other parts of the code;  is there a G4 bug that causes something to be
      // unitialized?
      visAtt->SetForceAuxEdgeVisible(forceAuxEdgeVisible);
      info.logical->SetVisAttributes(visAtt);

    }

    // Save the volume information in case someone else needs to access it by name.
    _helper.addVolInfo(info);

    return info;
  }

}
