//
// Free function to construct the stopping targets.
//
// $Id: constructStoppingTarget.cc,v 1.20 2013/08/30 17:00:14 genser Exp $
// $Author: genser $
// $Date: 2013/08/30 17:00:14 $
//
// Original author Peter Shanahan
//
// Notes:


// C++ includes
#include <iostream>
#include <string>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/constructStoppingTarget.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/checkForOverlaps.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"

using namespace std;

namespace mu2e {

  VolumeInfo constructStoppingTarget( VolumeInfo   const& parent,
                                      SimpleConfig const& config ){

    const bool forceAuxEdgeVisible  = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck       = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV              = true;

    int verbosity(config.getInt("target.verbosity",0));
    if ( verbosity > 1 ) std::cout<<"In constructStoppingTarget"<<std::endl;
    // Master geometry for the Target assembly
    GeomHandle<StoppingTarget> target;

    G4VSensitiveDetector* stSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::StoppingTarget());

    G4Helper    & _helper = *(art::ServiceHandle<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    TubsParams targetMotherParams(0., target->cylinderRadius(), target->cylinderLength()/2.);

    VolumeInfo targetInfo = nestTubs("StoppingTargetMother",
                                     targetMotherParams,
                                     findMaterialOrThrow(target->fillMaterial()),
                                     0,
                                     target->centerInMu2e() - parent.centerInMu2e(),
                                     parent,
                                     0,
                                     false/*visible*/,
                                     G4Colour::Black(),
                                     false/*solid*/,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

    // now create the individual targets

    G4VPhysicalVolume* pv;

    for (int itf=0; itf<target->nFoils(); itf++) {

        TargetFoil foil=target->foil(itf);

        VolumeInfo foilInfo;
        G4Material* foilMaterial = findMaterialOrThrow(foil.material());
        foilInfo.name = "Foil";

        foilInfo.solid = new G4Tubs(foilInfo.name
                                    ,foil.rIn()
                                    ,foil.rOut()
                                    ,foil.halfThickness()
                                    ,0.
                                    ,CLHEP::twopi
                                    );

        foilInfo.logical = new G4LogicalVolume( foilInfo.solid
                                                , foilMaterial
                                                , foilInfo.name
                                                );
        if(stSD) foilInfo.logical->SetSensitiveDetector(stSD);

        // rotation matrix...
        G4RotationMatrix* rot = 0; //... will have to wait

        G4ThreeVector foilOffset(foil.centerInMu2e() - targetInfo.centerInMu2e());
        if ( verbosity > 1 ) cout<<"foil "<<itf<<" centerInMu2e="<<foil.centerInMu2e()<<", offset="<<foilOffset<<endl;

        // G4 manages the lifetime of this object.
        pv = new G4PVPlacement( rot,
                           foilOffset,
                           foilInfo.logical,
                           "TargetFoil_",
                           targetInfo.logical,
                           0,
                           itf,
                           false);

        doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);

        if (!config.getBool("target.visible",true)) {
          foilInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
        } else {
          G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Magenta()));
          visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
          visAtt->SetForceSolid(config.getBool("target.solid",true));
          foilInfo.logical->SetVisAttributes(visAtt);
        }
      }// target foils

    return targetInfo;
  }

} // end namespace mu2e
