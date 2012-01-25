// Andrei Gaponenko, 2011

#include "Mu2eG4/inc/constructProtonBeamDump.hh"
#include "Mu2eG4/inc/constructExtMonFNAL.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"

#include <iostream>
#include <cmath>

#include "G4Color.hh"
#include "G4LogicalVolume.hh"
#include "G4Trap.hh"
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4Polycone.hh"
#include "G4ExtrudedSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4TwoVector.hh"
#include "G4Transform3D.hh"
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

#define AGDEBUG(stuff) std::cerr<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
//#define AGDEBUG(stuff)

namespace mu2e {

  //================================================================
  void constructCollimatorExtMonFNAL(const ProtonBeamDump::CollimatorExtMonFNAL& collimator, 
				     const VolumeInfo& parent,
				     const CLHEP::Hep3Vector& collimatorCenterInParent,
				     const SimpleConfig& config
				     )
  {
    GeomHandle<ProtonBeamDump> dump;

    MaterialFinder materialFinder(config);
    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV             = true;


    // Make sure the solids definining the collimator hole etc are sufficiently long to completely 
    // exit the concrete on both ends.

    const double boxdz = 0.5*collimator.horizontalLength();
    const double dr =  *std::max_element(collimator.alignmentPlugRadius().begin(), collimator.alignmentPlugRadius().end());
    const double cylHalfLength = std::sqrt(std::pow(boxdz, 2) +
					   std::pow(boxdz * tan(collimator.angleV()) + dr/cos(collimator.angleV()), 2) +
					   std::pow(boxdz * tan(collimator.angleH()) + dr/cos(collimator.angleH()), 2)
					   );

    double zPlane[] = {-cylHalfLength, -0.5*collimator.radiusTransitiondZ(), +0.5*collimator.radiusTransitiondZ(), +cylHalfLength };
    double rzero[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    
    G4Box* cutbox = reg.add(new G4Box(collimator.name()+"cutbox", 
				      dump->enclosureHalfSize()[0], 
				      dump->enclosureHalfSize()[1], 
				      0.5*collimator.horizontalLength())
			    );
    
    CLHEP::HepRotation *colrot = reg.add(CLHEP::HepRotation::IDENTITY);
    colrot->rotateX(-collimator.angleV());
    colrot->rotateY(-collimator.angleH());

    //----------------------------------------------------------------
    // Alignment hole
      
    double rAlignmentHole[] = {
      collimator.alignmentPlugRadius()[1] + collimator.alignmentHoleRClearance()[1],
      collimator.alignmentPlugRadius()[1] + collimator.alignmentHoleRClearance()[1],
      collimator.alignmentPlugRadius()[0] + collimator.alignmentHoleRClearance()[0],
      collimator.alignmentPlugRadius()[0] + collimator.alignmentHoleRClearance()[0]
    };

    G4Polycone *holeCylinder = reg.add(new G4Polycone(collimator.name()+"holecomponent", 0, 2*M_PI,
						      sizeof(zPlane)/sizeof(zPlane[0]),
						      zPlane,
						      rzero,
						      rAlignmentHole
						      )
				       );
      

    VolumeInfo alignmentHole(collimator.name()+"AlignmentHole",
			     collimatorCenterInParent,
			     parent.centerInWorld);
      
    alignmentHole.solid = new G4IntersectionSolid(alignmentHole.name,
						  holeCylinder,
						  cutbox,
						  G4Transform3D(*colrot, CLHEP::Hep3Vector(0,0,0))
						  );

    finishNesting(alignmentHole, 
		  materialFinder.get("hall.insideMaterialName"),
		  colrot,
		  alignmentHole.centerInParent,
		  parent.logical,
		  0, 
		  config.getBool("extMonFilter."+collimator.name()+".alignmentHole.visible"),
		  G4Colour::Cyan(),
		  config.getBool("extMonFilter."+collimator.name()+".alignmentHole.solid"),
		  forceAuxEdgeVisible,
		  placePV,
		  doSurfaceCheck
		  );

    //----------------------------------------------------------------
    // Alignment plug

    double rAlignmentPlug[] = {
      collimator.alignmentPlugRadius()[1],
      collimator.alignmentPlugRadius()[1],
      collimator.alignmentPlugRadius()[0],
      collimator.alignmentPlugRadius()[0]
    };

    G4Polycone *plugCylinder = reg.add(new G4Polycone(collimator.name()+"plugcomponent", 0, 2*M_PI,
						      sizeof(zPlane)/sizeof(zPlane[0]),
						      zPlane,
						      rzero,
						      rAlignmentPlug
						      )
				       );

    VolumeInfo alignmentPlug(collimator.name()+"AlignmentPlug",
			     CLHEP::Hep3Vector(0,0,0),
			     alignmentHole.centerInWorld);
      
    alignmentPlug.solid = new G4IntersectionSolid(alignmentPlug.name,
						  plugCylinder,
						  cutbox,
						  G4Transform3D(*colrot, CLHEP::Hep3Vector(0,0,0))
						  );
      
    finishNesting(alignmentPlug, 
		  materialFinder.get("protonBeamDump.material.shielding"),
		  0,
		  alignmentPlug.centerInParent,
		  alignmentHole.logical,
		  0, 
		  config.getBool("extMonFilter."+collimator.name()+".alignmentPlug.visible"),
		  G4Colour::Red(),
		  config.getBool("extMonFilter."+collimator.name()+".alignmentPlug.solid"),
		  forceAuxEdgeVisible,
		  placePV,
		  doSurfaceCheck
		  );

    //----------------------------------------------------------------
    // Collimator channel: the upstream half
    
    G4Box *upbox = reg.add(new G4Box(collimator.name()+"upbox", 
				     collimator.channelWidth()[0],
				     collimator.channelHeight()[0],
				     0.5*cylHalfLength
				     )
			   );
    
    VolumeInfo channelUp(collimator.name()+"ChannelUp",CLHEP::Hep3Vector(0,0,0.5*cylHalfLength), alignmentPlug.centerInWorld);
    
    channelUp.solid = new G4IntersectionSolid(channelUp.name,
					      upbox,
					      cutbox,
					      G4Transform3D(*colrot, CLHEP::Hep3Vector(0,0,-0.5*cylHalfLength))
					      );
    
    finishNesting(channelUp, 
		  materialFinder.get("hall.insideMaterialName"),
		  0,
		  channelUp.centerInParent,
		  alignmentPlug.logical,
		  0, 
		  config.getBool("extMonFilter."+collimator.name()+".channel.visible"),
		  G4Colour::Blue(),
		  config.getBool("extMonFilter."+collimator.name()+".channel.solid"),
		  forceAuxEdgeVisible,
		  placePV,
		  doSurfaceCheck
		  );
    
    //----------------------------------------------------------------
    // Collimator channel: the downstream half

    G4Box *collimatorDownBox = reg.add(new G4Box(collimator.name()+"dnbox", 
						 collimator.channelWidth()[1],
						 collimator.channelHeight()[1],
						 0.5*cylHalfLength)
				       );
      
    VolumeInfo collimatorDown(collimator.name()+"ChannelDn",CLHEP::Hep3Vector(0,0,-0.5*cylHalfLength), alignmentPlug.centerInWorld);
 
    collimatorDown.solid = new G4IntersectionSolid(collimatorDown.name,
						   collimatorDownBox,
						   cutbox,
						   G4Transform3D(*colrot, CLHEP::Hep3Vector(0,0,0.5*cylHalfLength))
						   );
      
    finishNesting(collimatorDown, 
		  materialFinder.get("hall.insideMaterialName"),
		  0,
		  collimatorDown.centerInParent,
		  alignmentPlug.logical,
		  0, 
		  config.getBool("extMonFilter."+collimator.name()+".channel.visible"),
		  G4Colour::Blue(),
		  config.getBool("extMonFilter."+collimator.name()+".channel.solid"),
		  forceAuxEdgeVisible,
		  placePV,
		  doSurfaceCheck
		  );

  }

  //================================================================
  void constructFilterMagnet(const ProtonBeamDump& dump, 
			     const VolumeInfo& parent,
			     const SimpleConfig& config
			     )
  {
    MaterialFinder materialFinder(config);

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV             = true;

    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

    //----------------------------------------------------------------
    CLHEP::HepRotation *magnetRotationInParent = reg.add(CLHEP::HepRotation::IDENTITY);
    magnetRotationInParent->rotateX(-dump.filterMagnetAngleV());

    const VolumeInfo magnetIron = nestBox("ExtMonFNALFilterMagnetIron",
					  dump.filterMagnet().outerHalfSize(),
					  materialFinder.get("extMonFilter.magnet.material"),
					  magnetRotationInParent,
					  dump.filterMagnetCenterInEnclosure() - dump.magnetPitCenterInEnclosure(),
					  parent, 0,
					  config.getBool("extMonFilter.magnet.iron.visible", true),
					  G4Colour::Magenta(),
					  config.getBool("extMonFilter.magnet.iron.solid", false),
					  forceAuxEdgeVisible,
					  placePV,
					  doSurfaceCheck
					  );

    std::vector<double> apertureHalfSize(3);
    apertureHalfSize[0] = 0.5*dump.filterMagnet().apertureWidth();
    apertureHalfSize[1] = 0.5*dump.filterMagnet().apertureHeight();
    apertureHalfSize[2] = dump.filterMagnet().outerHalfSize()[2];

    nestBox("ExtMonFNALFilterMagnetAperture",
	    apertureHalfSize,
	    materialFinder.get("hall.insideMaterialName"),
	    0,
	    CLHEP::Hep3Vector(0, 0, 0),
	    magnetIron.logical, 0,
	    config.getBool("extMonFilter.magnet.air.visible", true),
	    G4Colour::Blue(),
	    config.getBool("extMonFilter.magnet.air.solid", false),
	    forceAuxEdgeVisible,
	    placePV,
	    doSurfaceCheck
	    );
  }
    
  //================================================================

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
 
    const VolumeInfo magnetPit  = nestBox("ProtonBeamDumpMagnetPit",
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
    // The ExtMon filter channel

    constructCollimatorExtMonFNAL(dump->collimator1(),
				  logicalEnclosure,
				  dump->collimator1CenterInEnclosure(),
				  config);

    constructFilterMagnet(*dump, magnetPit, config);

    constructCollimatorExtMonFNAL(dump->collimator2(),
				  logicalEnclosure,
				  dump->collimator2CenterInEnclosure(),
				  config);

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
