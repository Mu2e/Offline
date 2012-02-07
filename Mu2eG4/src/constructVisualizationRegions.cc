// Andrei Gaponenko, 2011

#include "Mu2eG4/inc/constructVisualizationRegions.hh"

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include "G4Color.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Transform3D.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib/exception.h"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"

#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"

#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"


//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  //================================================================
  void constructVisualizationRegions(const VolumeInfo& worldVolume, const SimpleConfig& config)
  {
    GeomHandle<WorldG4> worldGeom;

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = false; // overlaps are OK
    const bool placePV             = true;

    std::vector<int> visible;
    config.getVectorInt("visregions.boxes.visible", visible, visible);
    if(!visible.empty()) {

      const int nboxes = visible.size();

      std::vector<int> solid; config.getVectorInt("visregions.boxes.solid", solid, nboxes);
      std::vector<std::string> material; config.getVectorString("visregions.boxes.material", material, nboxes);
      
      std::vector<double> red; config.getVectorDouble("visregions.boxes.color.red", red, nboxes);
      std::vector<double> green; config.getVectorDouble("visregions.boxes.color.green", green, nboxes);
      std::vector<double> blue; config.getVectorDouble("visregions.boxes.color.blue", blue, nboxes);


      std::vector<double> xmin; config.getVectorDouble("visregions.boxes.xmin", xmin, nboxes);
      std::vector<double> ymin; config.getVectorDouble("visregions.boxes.ymin", ymin, nboxes);
      std::vector<double> zmin; config.getVectorDouble("visregions.boxes.zmin", zmin, nboxes);

      std::vector<double> xmax; config.getVectorDouble("visregions.boxes.xmax", xmax, nboxes);
      std::vector<double> ymax; config.getVectorDouble("visregions.boxes.ymax", ymax, nboxes);
      std::vector<double> zmax; config.getVectorDouble("visregions.boxes.zmax", zmax, nboxes);

      for(int ibox = 0; ibox < nboxes; ++ibox) {
        
        CLHEP::Hep3Vector boxCenterInMu2e((xmax[ibox]+xmin[ibox])/2,
                                          (ymax[ibox]+ymin[ibox])/2,
                                          (zmax[ibox]+zmin[ibox])/2
                                          );

        std::vector<double> boxHalfSize(3);
        boxHalfSize[0] = (xmax[ibox]-xmin[ibox])/2;
        boxHalfSize[1] = (ymax[ibox]-ymin[ibox])/2;
        boxHalfSize[2] = (zmax[ibox]-zmin[ibox])/2;

        std::ostringstream boxname;
        boxname<<"VisualizationBox"<<std::setw(3)<<std::setfill('0')<<ibox;

        nestBox(boxname.str(),
                boxHalfSize,
                findMaterialOrThrow(material[ibox]),
                0,
                boxCenterInMu2e + worldGeom->mu2eOriginInWorld(),
                worldVolume, 0,
                visible[ibox],
                G4Colour(red[ibox], green[ibox], blue[ibox]),
                solid[ibox],
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
 	    );
      }

    }
//    //----------------------------------------------------------------
//    CLHEP::Hep3Vector enclosurePositionInDirt( beamDumpDirtRotation * (dump->enclosureCenterInMu2e() - beamDumpDirt.centerInMu2e()));
//
//    static CLHEP::HepRotation rotationInDirt(CLHEP::HepRotation::IDENTITY);
//    rotationInDirt.rotateZ(dump->coreRotY());
//    rotationInDirt.rotateX(+90*CLHEP::degree);
//
//    const VolumeInfo logicalEnclosure = nestBox("ProtonBeamDumpShielding",
// 						dump->enclosureHalfSize(), 
// 						materialFinder.get("protonBeamDump.material.shielding"),
//						&rotationInDirt,
// 						enclosurePositionInDirt,
// 						beamDumpDirt, 0,
//						config.getBool("protonBeamDump.logicalEnclosureVisible", true),
// 						G4Colour::Grey(),
//						config.getBool("protonBeamDump.logicalEnclosureSolid", false),
//						forceAuxEdgeVisible,
//						placePV,
//						doSurfaceCheck
// 						);
// 
//    nestBox("ProtonBeamDumpCore",
// 	    dump->coreHalfSize(),
// 	    materialFinder.get("protonBeamDump.material.core"),
// 	    0,
// 	    dump->coreCenterInEnclosure(),
// 	    logicalEnclosure, 0,
//	    config.getBool("protonBeamDump.coreVisible", true),
// 	    G4Colour::Blue(),
//	    config.getBool("protonBeamDump.coreSolid", true),
//	    forceAuxEdgeVisible,
//	    placePV,
//	    doSurfaceCheck
//
// 	    );
// 
//    nestBox("ProtonBeamDumpMouth",
// 	    dump->mouthHalfSize(),
// 	    materialFinder.get("protonBeamDump.material.air"),
// 	    0,
// 	    CLHEP::Hep3Vector(0.,
// 			      dump->coreCenterInEnclosure()[1],
// 			      dump->enclosureHalfSize()[2] - dump->mouthHalfSize()[2]),
// 	    logicalEnclosure, 0,
//	    config.getBool("protonBeamDump.mouthVisible", true), 
// 	    G4Colour::Cyan(), 
//	    config.getBool("protonBeamDump.mouthSolid", false), 
//	    forceAuxEdgeVisible,
//	    placePV,
//	    doSurfaceCheck
//
// 	    );
// 
//    nestBox("ProtonBeamNeutronCave",
// 	    dump->neutronCaveHalfSize(),
// 	    materialFinder.get("protonBeamDump.material.air"),
// 	    0,
// 	    CLHEP::Hep3Vector(0.,
// 			      dump->coreCenterInEnclosure()[1],
// 			      dump->enclosureHalfSize()[2] - 2*dump->mouthHalfSize()[2] - dump->neutronCaveHalfSize()[2]),
// 	    logicalEnclosure, 0,
//	    config.getBool("protonBeamDump.neutronCaveVisible", true), 
// 	    G4Colour::Cyan(),
//	    config.getBool("protonBeamDump.neutronCaveSolid", false), 
//	    forceAuxEdgeVisible,
//	    placePV,
//	    doSurfaceCheck
// 	    );
// 
//    const VolumeInfo magnetPit  = nestBox("ProtonBeamDumpMagnetPit",
//					  dump->magnetPitHalfSize(),
//					  materialFinder.get("protonBeamDump.material.air"),
//					  0,
//					  dump->magnetPitCenterInEnclosure(),
//					  logicalEnclosure, 0,
//					  config.getBool("protonBeamDump.magnetPitVisible", true), 
//					  G4Colour::Cyan(),
//					  config.getBool("protonBeamDump.magnetPitSolid", false),
//					  forceAuxEdgeVisible,
//					  placePV,
//					  doSurfaceCheck
//					  );
//
//    //----------------------------------------------------------------
//    // The ExtMon filter channel
//
//    constructCollimatorExtMonFNAL(dump->collimator1(),
//				  logicalEnclosure,
//				  dump->collimator1CenterInEnclosure(),
//				  config);
//
//    constructFilterMagnet(*dump, magnetPit, config);
//
//    constructCollimatorExtMonFNAL(dump->collimator2(),
//				  logicalEnclosure,
//				  dump->collimator2CenterInEnclosure(),
//				  config);
//
//    //----------------------------------------------------------------
//    // Create hall walls that are inside the formal hall box
//    
//    static CLHEP::HepRotation wallRotation(CLHEP::HepRotation::IDENTITY);
//    wallRotation.rotateX(-90*CLHEP::degree);
//    
//    // part on the +x side
//    if(true) {
//      std::vector<G4TwoVector> outline(building->concreteOuterOutline1());
//      outline.push_back(G4TwoVector(building->hallInsideXmax(), building->hallInsideZExtMonUCIWall() - building->hallWallThickness()));
//      outline.push_back(G4TwoVector(building->hallInsideXmax(), building->hallInsideZExtMonUCIWall()));
//      outline.push_back(G4TwoVector(dump->shieldingFaceXmax(), building->hallInsideZExtMonUCIWall()));
//
//      VolumeInfo wall("HallConcreteExtMonUCIWall",
//		      CLHEP::Hep3Vector(0, 0, 0),
//		      beamDumpDirt.centerInWorld);
//      
//      wall.solid = new G4ExtrudedSolid(wall.name, outline, 
//				       (building->hallInsideYmax()-building->hallInsideYmin())/2, 
//				       G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
//      
//      finishNesting(wall,
//		    materialFinder.get("hall.wallMaterialName"),
//		    0, //&wallRotation,
//		    wall.centerInParent,
//		    beamDumpDirt.logical,
//		    0, 
//		    config.getBool("hall.wallsVisible",true),
//		    G4Colour::Grey(),
//		    config.getBool("hall.wallsSolid",false),
//		    forceAuxEdgeVisible,
//		    placePV,
//		    doSurfaceCheck
//		    );
//    }
//
//    // part on the -x side
//    if(true) {
//      std::vector<G4TwoVector> outline(building->concreteOuterOutline3());
//      outline.push_back(G4TwoVector(dump->shieldingFaceXmin(), dump->shieldingFaceZatXmin()));
//      outline.push_back(G4TwoVector(building->hallInsideXmaxAtBeamDumpWall(), building->hallInsideZBeamDumpWall()));
//      outline.push_back(G4TwoVector(building->hallInsideXmin(), building->hallInsideZBeamDumpWall()));
//      outline.push_back(G4TwoVector(building->hallInsideXmin(), building->hallInsideZBeamDumpWall() - building->hallWallThickness()));
//
//      VolumeInfo wall("HallConcreteBeamDumpWall",
//		      CLHEP::Hep3Vector(0, 0, 0),
//		      beamDumpDirt.centerInWorld);
//      
//      wall.solid = new G4ExtrudedSolid(wall.name, outline, 
//				       (building->hallInsideYmax()-building->hallInsideYmin())/2, 
//				       G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
//      
//      finishNesting(wall,
//		    materialFinder.get("hall.wallMaterialName"),
//		    0,
//		    wall.centerInParent,
//		    beamDumpDirt.logical,
//		    0, 
//		    config.getBool("hall.wallsVisible",true),
//		    G4Colour::Grey(),
//		    config.getBool("hall.wallsSolid",false),
//		    forceAuxEdgeVisible,
//		    placePV,
//		    doSurfaceCheck
//		    );
//    }
//
//    //----------------------------------------------------------------
//    // Add some volumes for visualization purposes.
//    // 
//    // ROOT's OpenGL graphics refuses to display in a nice way the
//    // ProtonBeamDumpShielding volume.  On the other hand it turns out
//    // the same viewer shows existing VirtualDetectors nicely.  
//    // 
//    // Here we outline the dump shielding using very thin boxes (with
//    // the same thickness as the VirtualDetector dimension)
//    // 
//
//    const bool applyVisualizationKludge = config.getBool("protonBeamDump.applyROOTVisualizationKludge", false);
//    if(applyVisualizationKludge) {
//
//      const double kludgeHalfThickness = 0.01; // the thickness that works with the current root opengl 
//      int const nSurfaceCheckPoints = 100000; // for a more thorrow check due to the small thickness
//
//      const bool kludgeIsVisible      = true;
//      const bool kludgeIsSolid        = true; //false
//
//      G4Material* vacuumMaterial      = materialFinder.get("toyDS.insideMaterialName");
//
//      // The vertical side walls go inside the dump concrete
//      if(1) {
//
//	std::vector<double> hlen(3);
//	hlen[0] = kludgeHalfThickness;
//	hlen[1] = dump->enclosureHalfSize()[1];
//	hlen[2] = dump->enclosureHalfSize()[2];
//	
//	VolumeInfo pxInfo = nestBox("BeamDumpShieldingVisKludgePosX", 
//				    hlen,
//				    vacuumMaterial,
//				    0,
//				    CLHEP::Hep3Vector(+dump->enclosureHalfSize()[0] - kludgeHalfThickness, 0., 0.),
//				    logicalEnclosure,
//				    0,
//				    kludgeIsVisible,
//				    G4Color::Grey(),
//				    kludgeIsSolid,
//				    forceAuxEdgeVisible,
//				    true,
//				    false /*surface check*/
//				    );
//
//	// the volumes are very thin, a more thorough check is needed
//	doSurfaceCheck && pxInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
//
//	VolumeInfo nxInfo = nestBox("BeamDumpShieldingVisKludgeNegX", 
//				    hlen,
//				    vacuumMaterial,
//				    0,
//				    CLHEP::Hep3Vector(-dump->enclosureHalfSize()[0] + kludgeHalfThickness, 0., 0.),
//				    logicalEnclosure,
//				    0,
//				    kludgeIsVisible,
//				    G4Color::Grey(),
//				    kludgeIsSolid,
//				    forceAuxEdgeVisible,
//				    true,
//				    false /*surface check*/
//				    );
//	
//	// the volumes are very thin, a more thorough check is needed
//	doSurfaceCheck && nxInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
//      }
//
//      // The top "visualization plane kludge" can't be put inside because of the magnet pit volume
//      // So put both top and bottom planes outside
//      // 
//      // The top surface:
//      if(1) {
//
//	std::vector<double> hlen(3);
//	hlen[0] = dump->enclosureHalfSize()[0];
//	hlen[1] = kludgeHalfThickness;
//	hlen[2] = dump->enclosureHalfSize()[2];
//
//	CLHEP::Hep3Vector kludgeCenterInMu2e(dump->enclosureCenterInMu2e()[0],
//					     dump->enclosureCenterInMu2e()[1]+dump->enclosureHalfSize()[1],
//					     dump->enclosureCenterInMu2e()[2]
//					     );
//	
//	CLHEP::Hep3Vector kludgeCenterInDirt( beamDumpDirtRotation*(kludgeCenterInMu2e - beamDumpDirt.centerInMu2e()) );
//	
//	
//	CLHEP::Hep3Vector kludgeOffset(0, 0, kludgeHalfThickness);
//	
//	VolumeInfo pyInfo = nestBox("BeamDumpShieldingVisKludgePosY", 
//				    hlen,
//				    vacuumMaterial,
//				    &rotationInDirt,
//				    kludgeCenterInDirt + kludgeOffset,
//				    beamDumpDirt, // logicalEnclosure,
//				    0,
//				    kludgeIsVisible,
//				    G4Color::Grey(),
//				    kludgeIsSolid,
//				    forceAuxEdgeVisible,
//				    true,
//				    false /*surface check*/
//				    );
//
//	// the volumes are very thin, a more thorough check is needed
//	doSurfaceCheck && pyInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
//      }
//      
//      // The bottom surface:
//      if(1) {
//
//	std::vector<double> hlen(3);
//	hlen[0] = dump->enclosureHalfSize()[0];
//	hlen[1] = kludgeHalfThickness;
//	hlen[2] = dump->enclosureHalfSize()[2];
//
//	CLHEP::Hep3Vector kludgeCenterInMu2e(dump->enclosureCenterInMu2e()[0],
//					     dump->enclosureCenterInMu2e()[1]-dump->enclosureHalfSize()[1],
//					     dump->enclosureCenterInMu2e()[2]
//					     );
//	
//	CLHEP::Hep3Vector kludgeCenterInDirt( beamDumpDirtRotation*(kludgeCenterInMu2e - beamDumpDirt.centerInMu2e()) );
//	
//	
//	CLHEP::Hep3Vector kludgeOffset(0, 0, -kludgeHalfThickness);
//	
//	VolumeInfo nyInfo = nestBox("BeamDumpShieldingVisKludgeNegY", 
//				    hlen,
//				    vacuumMaterial,
//				    &rotationInDirt,
//				    kludgeCenterInDirt + kludgeOffset,
//				    beamDumpDirt, // logicalEnclosure,
//				    0,
//				    kludgeIsVisible,
//				    G4Color::Grey(),
//				    kludgeIsSolid,
//				    forceAuxEdgeVisible,
//				    true,
//				    false /*surface check*/
//				    );
//
//	// the volumes are very thin, a more thorough check is needed
//	doSurfaceCheck && nyInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
//      }
//    }
//
//    //----------------------------------------------------------------
//    if(art::ServiceHandle<GeometryService>()->hasElement<mu2e::ExtMonFNAL::ExtMon>()) {
//      constructExtMonFNAL(beamDumpDirt, beamDumpDirtRotation, &rotationInDirt, config);
//    }

  }


}
