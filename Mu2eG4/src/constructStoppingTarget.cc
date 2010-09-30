//
// Free function to construct the stopping targets.
//
// $Id: constructStoppingTarget.cc,v 1.5 2010/09/30 17:32:46 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/09/30 17:32:46 $
//
// Original author Peter Shanahan
//
// Notes:


// C++ includes
#include <iostream>
#include <string>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/constructStoppingTarget.hh"
#include "TargetGeom/inc/Target.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestTubs.hh"

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

    double rOut  = CLHEP::mm * target->cylinderRadius();
    double zHalf = CLHEP::mm * target->cylinderLength()/2.;

    // center in detector coords, assumed to be on axis
    double z0    = CLHEP::mm * target->cylinderCenter();

    VolumeInfo targetInfo;

    // Make the mother volume for the Target
    string targetName("StoppingTargetMother");
    std::cout<<"Looking for material "<<target->fillMaterial()<<std::endl;
    G4Material* fillMaterial = findMaterialOrThrow(target->fillMaterial());
    std::cout<<"Done Looking for material "<<target->fillMaterial()<<std::endl;
    G4ThreeVector targetOffset(0.,0.,(12000+z0-zOff));

    cout << "Target Offset: z0, zOff, z0-zOff: " 
         << z0 << " "
         << zOff << " "
         << z0-zOff << " "
         << endl;
   


    targetInfo.solid  = new G4Tubs( targetName,
                                    0., rOut, zHalf, 0., 2.*M_PI );
    
    targetInfo.logical = new G4LogicalVolume( targetInfo.solid, fillMaterial, targetName); 
    
    std::cout<<"targetOffset="<<targetOffset<<std::endl;
    std::cout<<"mother has "<<mother->GetNoDaughters()<<" daughters"<<std::endl;
    std::cout<<" they are:"<<std::endl;
    for (int id=0; id<mother->GetNoDaughters(); id++) cout<<id<<"="<<
      mother->GetDaughter(id)->GetName()<<" at "<<mother->GetDaughter(id)->GetTranslation()<<std::endl;
    targetInfo.physical =  new G4PVPlacement( 0, 
                                              targetOffset, 
                                              targetInfo.logical, 
                                              targetName, 
                                              mother, 
                                              0, 
                                              0,
                                              config.getBool("g4.doSurfaceCheck",false));

    // Visualization attributes of the the mother volume.
    {
      // i.e., none...
      G4VisAttributes* visAtt = new G4VisAttributes(false);
      targetInfo.logical->SetVisAttributes(visAtt);
    }

    // now create the individual targets

    for (unsigned int itf=0; itf<target->nFoils(); itf++)
      {

        TargetFoil foil=target->foil(itf);
        VolumeInfo foilInfo;
        G4Material* foilMaterial = findMaterialOrThrow( foil.material() );
        string foilName("Foil");

        foilInfo.solid = new G4Tubs(foilName
                                    ,foil.rIn()
                                    ,foil.rOut()
                                    ,foil.halfThickness()
                                    ,0.
                                    ,CLHEP::twopi*CLHEP::radian
                                    );

        foilInfo.logical = new G4LogicalVolume( foilInfo.solid
                                                , foilMaterial
                                                , foilName
                                                );

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

        // leak?
        G4VisAttributes* visAtt = new G4VisAttributes(true, G4Colour::Magenta() );
        visAtt->SetVisibility(config.getBool("target.visible",true));
        visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
        visAtt->SetForceSolid(config.getBool("target.solid",true));
        foilInfo.logical->SetVisAttributes(visAtt);
       
      }// target foils

    return targetInfo;
  }


} // end namespace mu2e

