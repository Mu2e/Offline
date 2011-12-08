// Andrei Gaponenko, 2011

#include "Mu2eG4/inc/constructProtonBeamDump.hh"

#include <iostream>

#include "G4Color.hh"
#include "G4LogicalVolume.hh"
#include "G4Trap.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4SDManager.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/ProtonBeamDump.hh"
#include "GeometryService/inc/Mu2eBuilding.hh"

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
 
    MaterialFinder materialFinder(config);
     
    const VolumeInfo logicalEnclosure = nestBox("ProtonBeamDumpShielding",
 						dump->enclosureHalfSize(), 
 						materialFinder.get("protonBeamDump.material.shielding"),
 						&dump->enclosureRotationInMu2e(), // assume the parent is not rotated
 						dump->enclosureCenterInMu2e() - parent.centerInMu2e(),
 						parent, 0, config.getBool("protonBeamDump.logicalEnclosureVisible"),
 						G4Colour::Grey(), false, true, true, true
 						);
 
    nestBox("ProtonBeamDumpCore",
 	    dump->coreHalfSize(),
 	    materialFinder.get("protonBeamDump.material.core"),
 	    0,
 	    dump->coreCenterInEnclosure(),
 	    logicalEnclosure, 0, true, 
 	    G4Colour::Blue(), true, true, true, true
 	    );
 
    nestBox("ProtonBeamDumpMouth",
 	    dump->mouthHalfSize(),
 	    materialFinder.get("protonBeamDump.material.air"),
 	    0,
 	    CLHEP::Hep3Vector(0.,
 			      dump->coreCenterInEnclosure()[1],
 			      dump->enclosureHalfSize()[2] - dump->mouthHalfSize()[2]),
 	    logicalEnclosure, 0, true, 
 	    G4Colour::Cyan(), false, true, true, true
 	    );
 
    nestBox("ProtonBeamNeutronCave",
 	    dump->neutronCaveHalfSize(),
 	    materialFinder.get("protonBeamDump.material.air"),
 	    0,
 	    CLHEP::Hep3Vector(0.,
 			      dump->coreCenterInEnclosure()[1],
 			      dump->enclosureHalfSize()[2] - 2*dump->mouthHalfSize()[2] - dump->neutronCaveHalfSize()[2]),
 	    logicalEnclosure, 0, true, 
 	    G4Colour::Cyan(), false, true, true, true
 	    );
 
    nestBox("ProtonBeamDumpMagnetPit",
 	    dump->magnetPitHalfSize(),
 	    materialFinder.get("protonBeamDump.material.air"),
 	    0,
 	    dump->magnetPitCenterInEnclosure(),
 	    logicalEnclosure, 0, true, 
 	    G4Colour::Cyan(), false, true, true, true
 	    );
 
//tmp:     //----------------------------------------------------------------
//tmp:     // Compute a trapezoid to connect the dump enclosure to the hall air
//tmp: 
//tmp:     std::cerr<<"AG: old hallAirZmin = "<<building->hallCenterInMu2e()[2] - building->hallInsideHalfLengths()[2]<<std::endl;
//tmp: 
//tmp:     const double hallAirZmin = -12*CLHEP::meter;// building->hallCenterInMu2e()[2] - building->hallInsideHalfLengths()[2];
//tmp: 
//tmp:     // distance from the enclosure center to the plane z=z_min(hall air)
//tmp:     // along the beam dump axis
//tmp:     const double lcenter = (hallAirZmin - dump->enclosureCenterInMu2e()[2])/std::abs(cos(dump->coreRotY()));
//tmp: 
//tmp:     AGDEBUG("hallAirZmin = "<<hallAirZmin<<", dz = "<<(hallAirZmin - dump->enclosureCenterInMu2e()[2])<<", lcenter = "<<lcenter);
//tmp: 
//tmp:     // the wedge sides are from the dump enclosure face to the z=z_min(hall air) plane,
//tmp:     // parallel to the dump axis.  Their lengths are
//tmp:     const double wdxLong  = lcenter - dump->enclosureHalfSize()[2] + dump->enclosureHalfSize()[0]*std::abs(tan(dump->coreRotY()));
//tmp:     const double wdxShort = lcenter - dump->enclosureHalfSize()[2] - dump->enclosureHalfSize()[0]*std::abs(tan(dump->coreRotY()));
//tmp: 
//tmp:     AGDEBUG("enclosureHalfSize[2] = "<<dump->enclosureHalfSize()[2]<<", coreRotY = "<<dump->coreRotY()<<", tan = "<<tan(dump->coreRotY()));
//tmp:     AGDEBUG("wdxLong = "<<wdxLong<<", wdxShort = "<<wdxShort);
//tmp: 
//tmp:     // other wedge dimensions.  Note the wedge constructor we use takes full sizes, not half sizes.
//tmp:     const double wdy = 2*dump->enclosureHalfSize()[0];
//tmp:     const double wdz = 2*dump->enclosureHalfSize()[1];
//tmp: 
//tmp:     AGDEBUG("wdz = "<<wdz<<", wdy = "<<wdy);
//tmp: 
//tmp:     static CLHEP::HepRotation wedgeRotation(CLHEP::HepRotation::IDENTITY);
//tmp:     wedgeRotation.rotateX(-90*CLHEP::degree).rotateZ(-90*CLHEP::degree + dump->coreRotY());
//tmp: 
//tmp:     CLHEP::Hep3Vector wedgeCenterInMu2e;
//tmp: 
//tmp:     wedgeCenterInMu2e[0] = dump->enclosureCenterInMu2e()[0]
//tmp:       + sin(dump->coreRotY())*(dump->enclosureHalfSize()[2] + 0.25 * (wdxShort+wdxLong));
//tmp: 
//tmp:     wedgeCenterInMu2e[1] = dump->enclosureCenterInMu2e()[1];
//tmp: 
//tmp:     wedgeCenterInMu2e[2] = dump->enclosureCenterInMu2e()[2] 
//tmp:       + cos(dump->coreRotY())*(dump->enclosureHalfSize()[2] + 0.25 * (wdxShort+wdxLong));
//tmp:     
//tmp:     AGDEBUG("wedgeCenterInMu2e = "<<wedgeCenterInMu2e);
//tmp: 
//tmp:     VolumeInfo wi;
//tmp:     wi.name = "PBDWedge";
//tmp:     wi.solid = new G4Trap(wi.name, wdz, wdy, wdxLong, wdxShort);
//tmp: 
//tmp:     finishNesting(wi, 
//tmp: 		  materialFinder.get("protonBeamDump.material.air"),
//tmp: 		  &wedgeRotation,
//tmp: 		  wedgeCenterInMu2e - parent.centerInMu2e(),
//tmp: 		  parent.logical,
//tmp: 		  0,
//tmp: 		  config.getBool("protonBeamDump.logicalEnclosureVisible"),
//tmp: 		  G4Colour::Cyan(),
//tmp: 		  false, true, true, true
//tmp: 		  );
//tmp: 
//tmp: //test:    //----------------------------------------------------------------
//tmp: //test:    // Create the "-z" concrete wall of the hall around the air wedge
//tmp: //test:
//tmp: //test:    const double hallNegZWallCenterZ = hallAirZmin - 0.5 * building->hallWallThickness();
//tmp: //test:    const double wedgeYmin = dump->enclosureCenterInMu2e()[1] - dump->enclosureHalfSize()[1];
//tmp: //test:
//tmp: //test:    const double bottomSlabHeight = wedgeYmin - 
//tmp: //test:      (building->hallCenterInMu2e()[1] - building->hallInsideHalfLengths()[1] - building->hallFloorThickness());
//tmp: //test:
//tmp: //test:    std::vector<double> bottomSlabHalfSize(3);
//tmp: //test:    bottomSlabHalfSize[0] = building->hallInsideHalfLengths()[0] + building->hallWallThickness();
//tmp: //test:    bottomSlabHalfSize[1] = 0.5*bottomSlabHeight;
//tmp: //test:    bottomSlabHalfSize[2] = 0.5*building->hallWallThickness();
//tmp: //test:
//tmp: //test:    nestBox("HallNegZWallBottomSlab",
//tmp: //test:	    bottomSlabHalfSize,
//tmp: //test:	    materialFinder.get("hall.wallMaterialName"),
//tmp: //test:	    0,
//tmp: //test:	    CLHEP::Hep3Vector(building->hallCenterInMu2e()[0],
//tmp: //test:			      wedgeYmin - 0.5*bottomSlabHeight ,
//tmp: //test:			      hallNegZWallCenterZ)
//tmp: //test:	    - parent.centerInMu2e(),
//tmp: //test:	    parent, 0, config.getBool("hall.visible"),
//tmp: //test:	    G4Colour::Red(),
//tmp: //test:	    config.getBool("hall.solid",false), true, true, true
//tmp: //test:	    );
//tmp: //test:    
//tmp: //test:    const double wedgeYmax = wedgeYmin + 2*dump->enclosureHalfSize()[1];
//tmp: //test:    const double topSlabHeight = 
//tmp: //test:      (building->hallCenterInMu2e()[1] + building->hallInsideHalfLengths()[1] + building->hallCeilingThickness())
//tmp: //test:      - wedgeYmax;
//tmp: //test:
//tmp: //test:    std::vector<double> topSlabHalfSize(3);
//tmp: //test:    topSlabHalfSize[0] = building->hallInsideHalfLengths()[0] + building->hallWallThickness();
//tmp: //test:    topSlabHalfSize[1] = 0.5*topSlabHeight;
//tmp: //test:    topSlabHalfSize[2] = 0.5*building->hallWallThickness();
//tmp: //test:
//tmp: //test:    nestBox("HallNegZWallTopSlab",
//tmp: //test:	    topSlabHalfSize,
//tmp: //test:	    materialFinder.get("hall.wallMaterialName"),
//tmp: //test:	    0,
//tmp: //test:	    CLHEP::Hep3Vector(building->hallCenterInMu2e()[0],
//tmp: //test:			      wedgeYmax + 0.5*topSlabHeight ,
//tmp: //test:			      hallNegZWallCenterZ)
//tmp: //test:	    - parent.centerInMu2e(),
//tmp: //test:	    parent, 0, config.getBool("hall.visible"),
//tmp: //test:	    G4Colour::Red(),
//tmp: //test:	    config.getBool("hall.solid",false), true, true, true
//tmp: //test:	    );
//tmp: //test:    
//tmp: //test:    //----------------
//tmp: //test:    // the opening in the wall for the beam dump
//tmp: //test:    const double xCenterOnInWall = dump->enclosureCenterInMu2e()[0] 
//tmp: //test:      + tan(dump->coreRotY()) * (hallAirZmin - dump->enclosureCenterInMu2e()[2]);
//tmp: //test:
//tmp: //test:
//tmp: //test:    const double openingHalfWidth = dump->enclosureHalfSize()[0]/cos(dump->coreRotY());
//tmp: //test:
//tmp: //test:    // the longer side of the "-x" slab
//tmp: //test:    const double snwdxLong = building->hallInsideHalfLengths()[0] + building->hallWallThickness()
//tmp: //test:      + (xCenterOnInWall - building->hallCenterInMu2e()[0]) - openingHalfWidth;
//tmp: //test:
//tmp: //test:    // the shorter side of the "+x" slab
//tmp: //test:    const double spwdxShort = building->hallInsideHalfLengths()[0] + building->hallWallThickness()
//tmp: //test:      - (xCenterOnInWall - building->hallCenterInMu2e()[0]) - openingHalfWidth;
//tmp: //test:
//tmp: //test:    // the difference between long and short sides
//tmp: //test:    const double sdx = building->hallWallThickness() * tan(dump->coreRotY());
//tmp: //test:
//tmp: //test:    // "-x" slab
//tmp: //test:    VolumeInfo snx;
//tmp: //test:    snx.name = "HallNegZWallNXSlab";
//tmp: //test:    snx.solid = new G4Trap(snx.name,
//tmp: //test:			   2*dump->enclosureHalfSize()[1],
//tmp: //test:			   building->hallWallThickness(),
//tmp: //test:			   snwdxLong, snwdxLong - sdx);
//tmp: //test:
//tmp: //test:    static CLHEP::HepRotation snxRotation(CLHEP::HepRotation::IDENTITY);
//tmp: //test:    snxRotation.rotateX(90*CLHEP::degree);
//tmp: //test:
//tmp: //test:    finishNesting(snx, 
//tmp: //test:		  materialFinder.get("hall.wallMaterialName"),
//tmp: //test:		  &snxRotation,
//tmp: //test:		  CLHEP::Hep3Vector(0.5*((xCenterOnInWall - openingHalfWidth - 0.5*sdx) + 
//tmp: //test:					 (building->hallCenterInMu2e()[0]
//tmp: //test:					  - building->hallInsideHalfLengths()[0] - building->hallWallThickness())
//tmp: //test:					 ),
//tmp: //test:				    dump->enclosureCenterInMu2e()[1],
//tmp: //test:				    hallAirZmin - 0.5*building->hallWallThickness())
//tmp: //test:		  - parent.centerInMu2e(),
//tmp: //test:		  parent.logical,
//tmp: //test:		  0,
//tmp: //test:		  config.getBool("protonBeamDump.logicalEnclosureVisible"),
//tmp: //test:		  G4Colour::Red(),
//tmp: //test:		  false, true, true, true
//tmp: //test:		  );
//tmp: //test:
//tmp: //test:    //----------------
//tmp: //test:    // "+x" slab
//tmp: //test:    VolumeInfo spx;
//tmp: //test:    spx.name = "HallNegZWallPXSlab";
//tmp: //test:    spx.solid = new G4Trap(spx.name,
//tmp: //test:			   2*dump->enclosureHalfSize()[1],
//tmp: //test:			   building->hallWallThickness(),
//tmp: //test:			   spwdxShort + sdx, spwdxShort);
//tmp: //test:
//tmp: //test:    static CLHEP::HepRotation spxRotation(CLHEP::HepRotation::IDENTITY);
//tmp: //test:    spxRotation.rotateY(180*CLHEP::degree).rotateX(90*CLHEP::degree);
//tmp: //test:
//tmp: //test:    finishNesting(spx, 
//tmp: //test:		  materialFinder.get("hall.wallMaterialName"),
//tmp: //test:		  &spxRotation,
//tmp: //test:		  CLHEP::Hep3Vector(0.5*((xCenterOnInWall + openingHalfWidth - 0.5*sdx) + 
//tmp: //test:					 (building->hallCenterInMu2e()[0]
//tmp: //test:					  + building->hallInsideHalfLengths()[0] + building->hallWallThickness())
//tmp: //test:					 ),
//tmp: //test:				    dump->enclosureCenterInMu2e()[1],
//tmp: //test:				    hallAirZmin - 0.5*building->hallWallThickness())
//tmp: //test:		  - parent.centerInMu2e(),
//tmp: //test:		  parent.logical,
//tmp: //test:		  0,
//tmp: //test:		  config.getBool("protonBeamDump.logicalEnclosureVisible"),
//tmp: //test:		  G4Colour::Red(),
//tmp: //test:		  false, true, true, true
//tmp: //test:		  );
//tmp: //test:
//tmp: //test:    //----------------------------------------------------------------
//tmp: //test:    // Visual references for debugging
//tmp: //test:    const bool addTestObjects  = false;
//tmp: //test:    if(addTestObjects) {
//tmp: //test:      VolumeInfo ti;
//tmp: //test:      ti.name = "TestSphere1";
//tmp: //test:      ti.solid = new G4Orb(ti.name, 0.1*CLHEP::meter);
//tmp: //test:      finishNesting(ti, 
//tmp: //test:		    materialFinder.get("protonBeamDump.material.air"),
//tmp: //test:		    0,
//tmp: //test:		    CLHEP::Hep3Vector(xCenterOnInWall,
//tmp: //test:				      dump->enclosureCenterInMu2e()[1],
//tmp: //test:				      hallAirZmin - 0.5*building->hallWallThickness())
//tmp: //test:		    - parent.centerInMu2e(),
//tmp: //test:		    parent.logical, 0, true, G4Colour::Red(), true, true, true, true
//tmp: //test:		    );
//tmp: //test:    }
//tmp: //test:    if(addTestObjects) {
//tmp: //test:      VolumeInfo ti;
//tmp: //test:      ti.name = "TestSphere2";
//tmp: //test:      ti.solid = new G4Orb(ti.name, 0.1*CLHEP::meter);
//tmp: //test:      finishNesting(ti, 
//tmp: //test:		    materialFinder.get("protonBeamDump.material.air"),
//tmp: //test:		    0,
//tmp: //test:		    CLHEP::Hep3Vector(0, 0, 0),
//tmp: //test:		    wi.logical, 0, true, G4Colour::Blue(), true, true, true, true
//tmp: //test:		    );
//tmp: //test:    }
//tmp: //test:    if(addTestObjects) {
//tmp: //test:      VolumeInfo ti;
//tmp: //test:      ti.name = "TestCyl1";
//tmp: //test:      ti.solid = new G4Tubs(ti.name, 0., 5*CLHEP::cm, openingHalfWidth, 0., CLHEP::twopi);
//tmp: //test:      static CLHEP::HepRotation r(CLHEP::HepRotation::IDENTITY);
//tmp: //test:      r.rotateY(90*CLHEP::degree);
//tmp: //test:      
//tmp: //test:      finishNesting(ti, 
//tmp: //test:		    materialFinder.get("protonBeamDump.material.air"),
//tmp: //test:		    &r,
//tmp: //test:		    CLHEP::Hep3Vector(xCenterOnInWall,
//tmp: //test:				      dump->enclosureCenterInMu2e()[1],
//tmp: //test:				      hallAirZmin)
//tmp: //test:		    - parent.centerInMu2e(),
//tmp: //test:		    parent.logical, 0, true, G4Colour::Red(), true, true, true, true
//tmp: //test:		    );
//tmp: //test:    }
//tmp: //test:    if(addTestObjects) {
//tmp: //test:      VolumeInfo ti;
//tmp: //test:      ti.name = "TestCylPBDE";
//tmp: //test:      ti.solid = new G4Tubs(ti.name, 0., 5*CLHEP::cm, lcenter, 0., CLHEP::twopi);
//tmp: //test:      finishNesting(ti, 
//tmp: //test:		    materialFinder.get("protonBeamDump.material.air"),
//tmp: //test:		    0,
//tmp: //test:		    CLHEP::Hep3Vector(0, 0, 0),
//tmp: //test:		    logicalEnclosure.logical, 0, true, G4Colour::Blue(), true, true, true, true
//tmp: //test:		    );
//tmp: //test:    }
//tmp:  
    //----------------------------------------------------------------
  }
}
