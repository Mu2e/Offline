//
// Free function to construct the stopping targets.
//
// $Id: constructStoppingTarget.cc,v 1.12 2011/05/18 14:21:44 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/18 14:21:44 $
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
//#include "Mu2eG4/inc/nestTubs.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/StoppingTargetSD.hh"
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

  VolumeInfo constructStoppingTarget( G4LogicalVolume* mother,
                                      double zOff,
                                      SimpleConfig const& config ){

    std::cout<<"In constructStoppingTarget"<<std::endl;
    // Master geometry for the Target assembly
    GeomHandle<Target> target;

    G4VSensitiveDetector* stSD = G4SDManager::GetSDMpointer()->
      FindSensitiveDetector(SensitiveDetectorName::StoppingTarget());

    G4Helper    & _helper = *(art::ServiceHandle<G4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    double rOut  = CLHEP::mm * target->cylinderRadius();
    double zHalf = CLHEP::mm * target->cylinderLength()/2.;

    // center in detector coords, assumed to be on axis
    double z0    = CLHEP::mm * target->cylinderCenter();

    VolumeInfo targetInfo;

    // Make the mother volume for the Target
    targetInfo.name = "StoppingTargetMother";
    std::cout<<"Looking for material "<<target->fillMaterial()<<std::endl;
    G4Material* fillMaterial = findMaterialOrThrow(target->fillMaterial());
    std::cout<<"Done Looking for material "<<target->fillMaterial()<<std::endl;
    G4ThreeVector targetOffset(0.,0.,(12000+z0-zOff));

    cout << "Target Offset: z0, zOff, z0-zOff: "
         << z0 << " "
         << zOff << " "
         << z0-zOff << " "
         << endl;



    targetInfo.solid  = new G4Tubs( targetInfo.name,
                                    0., rOut, zHalf, 0., 2.*M_PI );

    targetInfo.logical = new G4LogicalVolume( targetInfo.solid, fillMaterial, targetInfo.name);

    std::cout<<"targetOffset="<<targetOffset<<std::endl;
    std::cout<<"mother has "<<mother->GetNoDaughters()<<" daughters"<<std::endl;
    std::cout<<" they are:"<<std::endl;
    for (int id=0; id<mother->GetNoDaughters(); id++) cout<<id<<"="<<
      mother->GetDaughter(id)->GetName()<<" at "<<mother->GetDaughter(id)->GetTranslation()<<std::endl;
    targetInfo.physical =  new G4PVPlacement( 0,
                                              targetOffset,
                                              targetInfo.logical,
                                              targetInfo.name,
                                              mother,
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
                                    ,CLHEP::twopi*CLHEP::radian
                                    );

        foilInfo.logical = new G4LogicalVolume( foilInfo.solid
                                                , foilMaterial
                                                , foilInfo.name
                                                );
	foilInfo.logical->SetSensitiveDetector(stSD);

        // rotation matrix...
        G4RotationMatrix* rot = 0; //... will have to wait

        G4ThreeVector foilOffset(foil.center()-G4ThreeVector(0.,0.,z0));
        cout<<"foil "<<itf<<" center="<<foil.center()<<", offset="<<foilOffset<<endl;

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

