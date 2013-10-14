//
// Free function to construct the stopping targets.
//
// $Id: constructStoppingTarget.cc,v 1.21 2013/10/14 23:57:32 roehrken Exp $
// $Author: roehrken $
// $Date: 2013/10/14 23:57:32 $
//
// Original author Peter Shanahan
//
// Notes:


// C++ includes
#include <iostream>
#include <string>

//boost includes
#include <boost/lexical_cast.hpp>

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

    // now create the individual target foils
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
//                           ("TargetFoil_" + boost::lexical_cast<std::string>(itf)).c_str(), // problems with Analyses/StoppingTarget00_module.cc. This module detects the stops and writes out the coordinates to muonPointFile. Module requires stopping target volumes to have the name "TargetFoil_". So numbering is currently possible here.
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

    for (int itf=0; itf<target->nSupportStructures(); itf++) {
	TargetFoilSupportStructure supportStructure=target->supportStructure(itf);
        //TargetFoil foil=target->foil(itf);

        //VolumeInfo foilInfo;
        VolumeInfo supportStructureInfo;

        G4Material* supportStructureMaterial = findMaterialOrThrow(supportStructure.material());
        supportStructureInfo.name = "FoilSupportStructure";

        supportStructureInfo.solid = new G4Tubs(supportStructureInfo.name
                                    ,0
                                    ,supportStructure.radius()
                                    ,supportStructure.length()/2. // G4Tubs useses half-lengths as this parameter
                                    ,0.
                                    ,CLHEP::twopi
                                    );

        supportStructureInfo.logical = new G4LogicalVolume( supportStructureInfo.solid
                                                , supportStructureMaterial
                                                , supportStructureInfo.name
                                                );
        if(stSD) supportStructureInfo.logical->SetSensitiveDetector(stSD);

	if ( verbosity > 1 ) std::cout << "supportStructure.support_id() = " << supportStructure.support_id() << "    target->nSupportStructures() = " << target->nSupportStructures() << "     target->nFoils() = " << target->nFoils() << "     supportStructure.length() = " << supportStructure.length() << std::endl;

        // rotation matrices to rotate the orientation of the supporting wires. First rotate into xy-plane by 90deg rotation around y-axis, then rotate within xy-plane by appropiate rotation around z-axis
        CLHEP::HepRotationY secRy(-M_PI/2.);
        CLHEP::HepRotationZ secRz( -supportStructure.support_id() * 360.*deg / (target->nSupportStructures()/target->nFoils()) - 90.*deg - supportStructure.angleOffset()*deg);
        G4RotationMatrix* supportStructure_rotMatrix = reg.add(G4RotationMatrix(secRy*secRz));

	if ( verbosity > 1 ) std::cout << "supportStructure_rotMatrix = " << *supportStructure_rotMatrix << std::endl;

        // vector where to place to support tube
        G4ThreeVector supportStructureOffset(supportStructure.centerInMu2e() - targetInfo.centerInMu2e()); // first find target center

        if ( verbosity > 1 ) cout<<"FoilSupportStructure "<<itf<<" centerInMu2e="<<supportStructure.centerInMu2e()<<", offset="<<supportStructureOffset<<endl;

        G4ThreeVector vector_supportStructure_Orientation( (supportStructure.length()/2.+supportStructure.foil_outer_radius()) * std::cos(supportStructure.support_id() * 360.*deg / (target->nSupportStructures()/target->nFoils()) + 90.*deg + supportStructure.angleOffset()*deg), (supportStructure.length()/2.+supportStructure.foil_outer_radius()) * std::sin(supportStructure.support_id() * 360.*deg / (target->nSupportStructures()/target->nFoils()) + 90.*deg + supportStructure.angleOffset()*deg), 0);
	if ( verbosity > 1 ) std::cout << "vector_supportStructure_Orientation = " << vector_supportStructure_Orientation << std::endl;

	supportStructureOffset += vector_supportStructure_Orientation; // second add vector to support wire tube center
	if ( verbosity > 1 ) std::cout << "supportStructureOffset += vector_supportStructure_Orientation = " << supportStructureOffset << std::endl;

        pv = new G4PVPlacement( supportStructure_rotMatrix,
                           supportStructureOffset,
                           supportStructureInfo.logical,
                           ("TargetSupportStructure_" + boost::lexical_cast<std::string>(itf)).c_str(),
                           targetInfo.logical,
                           0,
                           itf,
                           false);

        doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);

        if (!config.getBool("target.visible",true)) {
          supportStructureInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
        } else {
          G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Blue()));
          visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
          visAtt->SetForceSolid(config.getBool("target.solid",true));
          supportStructureInfo.logical->SetVisAttributes(visAtt);
        }
      }// target foils support structures

    return targetInfo;
  }

} // end namespace mu2e
