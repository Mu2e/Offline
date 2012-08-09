//
// Free function to construct the stopping targets.
//
// $Id: constructStoppingTarget.cc,v 1.16 2012/08/09 22:22:53 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/08/09 22:22:53 $
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
#include "TargetGeom/inc/Target.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

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

  VolumeInfo constructStoppingTarget( VolumeInfo   const& mother,
                                      SimpleConfig const& config ){

    int verbosity(config.getInt("target.verbosity",0));
    if ( verbosity > 1 ) std::cout<<"In constructStoppingTarget"<<std::endl;
    // Master geometry for the Target assembly
    GeomHandle<Target> target;

    G4VSensitiveDetector* stSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::StoppingTarget());

    G4Helper    & _helper = *(art::ServiceHandle<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    double rOut  = target->cylinderRadius();
    double zHalf = target->cylinderLength()/2.;

    // center in detector coords, assumed to be on axis
    double z0    = target->cylinderCenter();

    VolumeInfo targetInfo;

    // Make the mother volume for the Target
    targetInfo.name = "StoppingTargetMother";
    if ( verbosity > 1) std::cout<<"Looking for material "<<target->fillMaterial()<<std::endl;
    G4Material* fillMaterial = findMaterialOrThrow(target->fillMaterial());
    if ( verbosity > 1) std::cout<<"Done Looking for material "<<target->fillMaterial()<<std::endl;

    G4ThreeVector targetOffset(0.,0.,12000+z0- mother.centerInMu2e().z()); 

    if ( verbosity > 1 ) std::cout<<"targetOffset="<<targetOffset<<std::endl;

    targetInfo.solid  = new G4Tubs( targetInfo.name,
                                    0., rOut, zHalf, 0., 2.*M_PI );

    targetInfo.logical = new G4LogicalVolume( targetInfo.solid, fillMaterial, targetInfo.name);

    int nnd = mother.logical->GetNoDaughters();
    if ( verbosity > 1 ) {
      std::cout<<"mother has "<<nnd<<" daughters"<<std::endl;
      // well the mother has no daughters yet... but let's leave most of the code here for now
      if (nnd > 0 ) {
        std::cout<<" they are:"<<std::endl;
        for (int id=0; id<nnd; id++) cout<<id<<"="<<
          mother.logical->GetDaughter(id)->GetName()<<
          " at "<<mother.logical->GetDaughter(id)->GetTranslation()<<std::endl;
      }
    }
    targetInfo.physical =  new G4PVPlacement( 0,
                                              targetOffset,
                                              targetInfo.logical,
                                              targetInfo.name,
                                              mother.logical,
                                              0,
                                              0,
                                              config.getBool("g4.doSurfaceCheck",false));

    // Visualization attributes of the the mother volume.
    {
      // i.e., none...
      targetInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
    }

    // now create the individual targets

    for (int itf=0; itf<target->nFoils(); itf++)
      {

        TargetFoil foil=target->foil(itf);
        VolumeInfo foilInfo;
        G4Material* foilMaterial = findMaterialOrThrow( foil.material() );
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
        foilInfo.logical->SetSensitiveDetector(stSD);

        // rotation matrix...
        G4RotationMatrix* rot = 0; //... will have to wait

        G4ThreeVector foilOffset(foil.center()-G4ThreeVector(0.,0.,z0));
        if ( verbosity > 1 ) cout<<"foil "<<itf<<" center="<<foil.center()<<", offset="<<foilOffset<<endl;

        // G4 manages the lifetime of this object.
        new G4PVPlacement( rot,
                           foilOffset,
                           foilInfo.logical,
                           "TargetFoil_",
                           targetInfo.logical,
                           0,
                           itf,
                           config.getBool("g4.doSurfaceCheck",false));

        if (!config.getBool("target.visible",true)) {
          foilInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
        } else {
          G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Magenta()));
          visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
          visAtt->SetForceSolid(config.getBool("target.solid",true));
          foilInfo.logical->SetVisAttributes(visAtt);
        }
      }// target foils

    // Save the volume information in case someone else needs to access it by name.
    _helper.addVolInfo(targetInfo);

    return targetInfo;
  }


} // end namespace mu2e

