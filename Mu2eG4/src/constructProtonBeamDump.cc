// Andrei Gaponenko, 2011

#include "Mu2eG4/inc/constructProtonBeamDump.hh"
#include "Mu2eG4/inc/constructExtMonFNAL.hh"

#include <iostream>

#include "G4Color.hh"
#include "G4LogicalVolume.hh"
#include "G4Trap.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4ExtrudedSolid.hh"
#include "G4TwoVector.hh"
#include "G4SDManager.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib/exception.h"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/ProtonBeamDump.hh"
#include "GeometryService/inc/Mu2eBuilding.hh"
#include "GeometryService/inc/WorldG4.hh"

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"

#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

#include "Mu2eG4/inc/finishNesting.hh"
#include "G4Helper/inc/VolumeInfo.hh"

//#define AGDEBUG(stuff) std::cerr<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {
  void constructProtonBeamDump(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<ProtonBeamDump> dump;
    GeomHandle<Mu2eBuilding> building;
    GeomHandle<WorldG4> world;
 
    MaterialFinder materialFinder(config);

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV             = true;

    //----------------------------------------------------------------
    // Re-fill a part of the formal "HallAir" with dirt.
    // Use an extruded solid to have a properly angled facet 
    // at which the beam dump can be placed.

    // points need to be in the clock-wise order
    std::vector<G4TwoVector> beamDumpDirtOutiline;

    beamDumpDirtOutiline.push_back(G4TwoVector(building->hallInsideXmaxAtBeamDumpWall(), building->hallInsideZBeamDumpWall()));

    beamDumpDirtOutiline.push_back(G4TwoVector(building->hallInsideXmin(), building->hallInsideZBeamDumpWall()));

    beamDumpDirtOutiline.push_back(G4TwoVector(building->hallInsideXmin(), world->hallFormalZminInMu2e()));

    beamDumpDirtOutiline.push_back(G4TwoVector(building->hallInsideXmax(), world->hallFormalZminInMu2e()));

    beamDumpDirtOutiline.push_back(G4TwoVector(building->hallInsideXmax(), building->hallInsideZExtMonUCIWall()));

    beamDumpDirtOutiline.push_back(G4TwoVector(dump->shieldingFaceXmax(), building->hallInsideZExtMonUCIWall()));

    beamDumpDirtOutiline.push_back(G4TwoVector(dump->shieldingFaceXmax(), dump->shieldingFaceZatXmax()));

    if(dump->shieldingFaceZatXmin() > building->hallInsideZBeamDumpWall()) {
      throw cet::exception("GEOM")<<"constructProtonBeamDump(): hallInsideZBeamDumpWall conflicts with the proton dump enclosure\n";
    }

    // Don't add the last point if it coincides with the first
    if(dump->shieldingFaceZatXmin() < building->hallInsideZBeamDumpWall()) {
      beamDumpDirtOutiline.push_back(G4TwoVector(dump->shieldingFaceXmin(), dump->shieldingFaceZatXmin()));
    }

    static CLHEP::HepRotation beamDumpDirtRotation(CLHEP::HepRotation::IDENTITY);
    beamDumpDirtRotation.rotateX(-90*CLHEP::degree);

    VolumeInfo beamDumpDirt("ProtonBeamDumpDirt",
			    CLHEP::Hep3Vector(0, (building->hallInsideYmin() + building->hallInsideYmax())/2, 0)
			    - parent.centerInMu2e(),
			    parent.centerInWorld);
    
    beamDumpDirt.solid = new G4ExtrudedSolid(beamDumpDirt.name, beamDumpDirtOutiline, 
					     (building->hallInsideYmax() - building->hallInsideYmin())/2, 
					     G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

    finishNesting(beamDumpDirt, 
		  materialFinder.get("dirt.overburdenMaterialName"),
		  &beamDumpDirtRotation,
		  beamDumpDirt.centerInParent,
		  parent.logical,
		  0, 
		  config.getBool("protonBeamDump.dirtVisible", false),
		  G4Colour::Magenta(),
		  config.getBool("protonBeamDump.dirtSolid", false),
		  forceAuxEdgeVisible,
		  placePV,
		  doSurfaceCheck
		  );

    //----------------------------------------------------------------
    CLHEP::Hep3Vector enclosurePositionInDirt( beamDumpDirtRotation * (dump->enclosureCenterInMu2e() - beamDumpDirt.centerInMu2e()));

    static CLHEP::HepRotation rotationInDirt(CLHEP::HepRotation::IDENTITY);
    rotationInDirt.rotateZ(dump->coreRotY());
    rotationInDirt.rotateX(+90*CLHEP::degree);

    const VolumeInfo logicalEnclosure = nestBox("ProtonBeamDumpShielding",
 						dump->enclosureHalfSize(), 
 						materialFinder.get("protonBeamDump.material.shielding"),
						&rotationInDirt,
 						enclosurePositionInDirt,
 						beamDumpDirt, 0,
						config.getBool("protonBeamDump.logicalEnclosureVisible", true),
 						G4Colour::Grey(),
						config.getBool("protonBeamDump.logicalEnclosureSolid", false),
						forceAuxEdgeVisible,
						placePV,
						doSurfaceCheck
 						);
 
    nestBox("ProtonBeamDumpCore",
 	    dump->coreHalfSize(),
 	    materialFinder.get("protonBeamDump.material.core"),
 	    0,
 	    dump->coreCenterInEnclosure(),
 	    logicalEnclosure, 0,
	    config.getBool("protonBeamDump.coreVisible", true),
 	    G4Colour::Blue(),
	    config.getBool("protonBeamDump.coreSolid", true),
	    forceAuxEdgeVisible,
	    placePV,
	    doSurfaceCheck

 	    );
 
    nestBox("ProtonBeamDumpMouth",
 	    dump->mouthHalfSize(),
 	    materialFinder.get("protonBeamDump.material.air"),
 	    0,
 	    CLHEP::Hep3Vector(0.,
 			      dump->coreCenterInEnclosure()[1],
 			      dump->enclosureHalfSize()[2] - dump->mouthHalfSize()[2]),
 	    logicalEnclosure, 0,
	    config.getBool("protonBeamDump.mouthVisible", true), 
 	    G4Colour::Cyan(), 
	    config.getBool("protonBeamDump.mouthSolid", false), 
	    forceAuxEdgeVisible,
	    placePV,
	    doSurfaceCheck

 	    );
 
    nestBox("ProtonBeamNeutronCave",
 	    dump->neutronCaveHalfSize(),
 	    materialFinder.get("protonBeamDump.material.air"),
 	    0,
 	    CLHEP::Hep3Vector(0.,
 			      dump->coreCenterInEnclosure()[1],
 			      dump->enclosureHalfSize()[2] - 2*dump->mouthHalfSize()[2] - dump->neutronCaveHalfSize()[2]),
 	    logicalEnclosure, 0,
	    config.getBool("protonBeamDump.neutronCaveVisible", true), 
 	    G4Colour::Cyan(),
	    config.getBool("protonBeamDump.neutronCaveSolid", false), 
	    forceAuxEdgeVisible,
	    placePV,
	    doSurfaceCheck
 	    );
 
    nestBox("ProtonBeamDumpMagnetPit",
 	    dump->magnetPitHalfSize(),
 	    materialFinder.get("protonBeamDump.material.air"),
 	    0,
 	    dump->magnetPitCenterInEnclosure(),
 	    logicalEnclosure, 0,
	    config.getBool("protonBeamDump.magnetPitVisible", true), 
 	    G4Colour::Cyan(),
	    config.getBool("protonBeamDump.magnetPitSolid", false),
	    forceAuxEdgeVisible,
	    placePV,
	    doSurfaceCheck
 	    );

    //----------------------------------------------------------------
    // Create hall walls that are inside the formal hall box
    
    static CLHEP::HepRotation wallRotation(CLHEP::HepRotation::IDENTITY);
    wallRotation.rotateX(-90*CLHEP::degree);
    
    // part on the +x side
    if(true) {
      std::vector<G4TwoVector> outline(building->concreteOuterOutline1());
      outline.push_back(G4TwoVector(building->hallInsideXmax(), building->hallInsideZExtMonUCIWall() - building->hallWallThickness()));
      outline.push_back(G4TwoVector(building->hallInsideXmax(), building->hallInsideZExtMonUCIWall()));
      outline.push_back(G4TwoVector(dump->shieldingFaceXmax(), building->hallInsideZExtMonUCIWall()));

      VolumeInfo wall("HallConcreteExtMonUCIWall",
		      CLHEP::Hep3Vector(0, 0, 0),
		      beamDumpDirt.centerInWorld);
      
      wall.solid = new G4ExtrudedSolid(wall.name, outline, 
				       (building->hallInsideYmax()-building->hallInsideYmin())/2, 
				       G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
      
      finishNesting(wall,
		    materialFinder.get("hall.wallMaterialName"),
		    0, //&wallRotation,
		    wall.centerInParent,
		    beamDumpDirt.logical,
		    0, 
		    config.getBool("hall.wallsVisible",true),
		    G4Colour::Grey(),
		    config.getBool("hall.wallsSolid",false),
		    forceAuxEdgeVisible,
		    placePV,
		    doSurfaceCheck
		    );
    }

    // part on the -x side
    if(true) {
      std::vector<G4TwoVector> outline(building->concreteOuterOutline3());
      outline.push_back(G4TwoVector(dump->shieldingFaceXmin(), dump->shieldingFaceZatXmin()));
      outline.push_back(G4TwoVector(building->hallInsideXmaxAtBeamDumpWall(), building->hallInsideZBeamDumpWall()));
      outline.push_back(G4TwoVector(building->hallInsideXmin(), building->hallInsideZBeamDumpWall()));
      outline.push_back(G4TwoVector(building->hallInsideXmin(), building->hallInsideZBeamDumpWall() - building->hallWallThickness()));

      VolumeInfo wall("HallConcreteBeamDumpWall",
		      CLHEP::Hep3Vector(0, 0, 0),
		      beamDumpDirt.centerInWorld);
      
      wall.solid = new G4ExtrudedSolid(wall.name, outline, 
				       (building->hallInsideYmax()-building->hallInsideYmin())/2, 
				       G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
      
      finishNesting(wall,
		    materialFinder.get("hall.wallMaterialName"),
		    0,
		    wall.centerInParent,
		    beamDumpDirt.logical,
		    0, 
		    config.getBool("hall.wallsVisible",true),
		    G4Colour::Grey(),
		    config.getBool("hall.wallsSolid",false),
		    forceAuxEdgeVisible,
		    placePV,
		    doSurfaceCheck
		    );
    }

    //----------------------------------------------------------------
    if(art::ServiceHandle<GeometryService>()->hasElement<mu2e::ExtMonFNAL::ExtMon>()) {
      constructExtMonFNAL(beamDumpDirt, beamDumpDirtRotation, &rotationInDirt, config);
    }

  }
}
