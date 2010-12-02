//
// Free function to be used by the nest... functions
// 
// $Id: finishNesting.cc,v 1.1 2010/12/02 17:44:32 genser Exp $
// $Author: genser $ 
// $Date: 2010/12/02 17:44:32 $
//
// Original author KLG
//

// Mu2e includes
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/VolumeInfo.hh"
#include "Mu2eG4/inc/G4Helper.hh"

// G4 includes
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"


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
                     ) {

    G4Helper    & _helper = *(edm::Service<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    info.logical  = new G4LogicalVolume( info.solid, material, info.name); 

    info.physical  =  placePV ?
      info.physical = new G4PVPlacement( rot,
                                         offset, 
                                         info.logical, 
                                         info.name, 
                                         parent, 
                                         0, 
                                         copyNo, 
                                         doSurfCheck)
      :
      0;


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

    return;

  }

}
