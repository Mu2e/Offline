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

    //----------------------------------------------------------------
    // Compute a trapezoid to connect the dump enclosure to the hall air

    const double hallAirZmin = building->hallCenterInMu2e()[2] - building->hallInsideHalfLengths()[2];

    // distance from the enclosure center to the plane z=z_min(hall air)
    // along the beam dump axis
    const double lcenter = (hallAirZmin - dump->enclosureCenterInMu2e()[2])/std::abs(cos(dump->coreRotY()));

    AGDEBUG("hallAirZmin = "<<hallAirZmin<<", dz = "<<(hallAirZmin - dump->enclosureCenterInMu2e()[2])<<", lcenter = "<<lcenter);

    // the wedge sides are from the dump enclosure face to the z=z_min(hall air) plane,
    // parallel to the dump axis.  Their lengths are
    const double wdxLong  = lcenter - dump->enclosureHalfSize()[2] + dump->enclosureHalfSize()[0]*std::abs(tan(dump->coreRotY()));
    const double wdxShort = lcenter - dump->enclosureHalfSize()[2] - dump->enclosureHalfSize()[0]*std::abs(tan(dump->coreRotY()));

    AGDEBUG("enclosureHalfSize[2] = "<<dump->enclosureHalfSize()[2]<<", coreRotY = "<<dump->coreRotY()<<", tan = "<<tan(dump->coreRotY()));
    AGDEBUG("wdxLong = "<<wdxLong<<", wdxShort = "<<wdxShort);

    // other wedge dimensions.  Note the wedge constructor we use takes full sizes, not half sizes.
    const double wdy = 2*dump->enclosureHalfSize()[0];
    const double wdz = 2*dump->enclosureHalfSize()[1];

    AGDEBUG("wdz = "<<wdz<<", wdy = "<<wdy);

    static CLHEP::HepRotation wedgeRotation(CLHEP::HepRotation::IDENTITY);
    wedgeRotation.rotateX(-90*CLHEP::degree).rotateZ(-90*CLHEP::degree + dump->coreRotY());

    CLHEP::Hep3Vector wedgeCenterInMu2e;

    wedgeCenterInMu2e[0] = dump->enclosureCenterInMu2e()[0]
      + sin(dump->coreRotY())*(dump->enclosureHalfSize()[2] + 0.25 * (wdxShort+wdxLong));

    wedgeCenterInMu2e[1] = dump->enclosureCenterInMu2e()[1];

    wedgeCenterInMu2e[2] = dump->enclosureCenterInMu2e()[2] 
      + cos(dump->coreRotY())*(dump->enclosureHalfSize()[2] + 0.25 * (wdxShort+wdxLong));
    
    AGDEBUG("wedgeCenterInMu2e = "<<wedgeCenterInMu2e);

    VolumeInfo wi;
    wi.name = "PBDWedge";
    wi.solid = new G4Trap(wi.name, wdz, wdy, wdxLong, wdxShort);

    finishNesting(wi, 
		  materialFinder.get("protonBeamDump.material.air"),
		  &wedgeRotation,
		  wedgeCenterInMu2e - parent.centerInMu2e(),
		  parent.logical,
		  0,
		  config.getBool("protonBeamDump.logicalEnclosureVisible"),
		  G4Colour::Cyan(),
		  false, true, true, true
		  );

    //----------------------------------------------------------------
    // Create the "-z" concrete wall of the hall around the air wedge

    const double hallNegZWallCenterZ = hallAirZmin - 0.5 * building->hallWallThickness();
    const double wedgeYmin = dump->enclosureCenterInMu2e()[1] - dump->enclosureHalfSize()[1];

    const double bottomSlabHeight = wedgeYmin - 
      (building->hallCenterInMu2e()[1] - building->hallInsideHalfLengths()[1] - building->hallFloorThickness());

    std::vector<double> bottomSlabHalfSize(3);
    bottomSlabHalfSize[0] = building->hallInsideHalfLengths()[0] + building->hallWallThickness();
    bottomSlabHalfSize[1] = 0.5*bottomSlabHeight;
    bottomSlabHalfSize[2] = 0.5*building->hallWallThickness();

    nestBox("HallNegZWallBottomSlab",
	    bottomSlabHalfSize,
	    materialFinder.get("hall.wallMaterialName"),
	    0,
	    CLHEP::Hep3Vector(building->hallCenterInMu2e()[0],
			      wedgeYmin - 0.5*bottomSlabHeight ,
			      hallNegZWallCenterZ)
	    - parent.centerInMu2e(),
	    parent, 0, config.getBool("hall.visible"),
	    G4Colour::Red(),
	    config.getBool("hall.solid",false), true, true, true
	    );
    
    const double wedgeYmax = wedgeYmin + 2*dump->enclosureHalfSize()[1];
    const double topSlabHeight = 
      (building->hallCenterInMu2e()[1] + building->hallInsideHalfLengths()[1] + building->hallCeilingThickness())
      - wedgeYmax;

    std::vector<double> topSlabHalfSize(3);
    topSlabHalfSize[0] = building->hallInsideHalfLengths()[0] + building->hallWallThickness();
    topSlabHalfSize[1] = 0.5*topSlabHeight;
    topSlabHalfSize[2] = 0.5*building->hallWallThickness();

    nestBox("HallNegZWallTopSlab",
	    topSlabHalfSize,
	    materialFinder.get("hall.wallMaterialName"),
	    0,
	    CLHEP::Hep3Vector(building->hallCenterInMu2e()[0],
			      wedgeYmax + 0.5*topSlabHeight ,
			      hallNegZWallCenterZ)
	    - parent.centerInMu2e(),
	    parent, 0, config.getBool("hall.visible"),
	    G4Colour::Red(),
	    config.getBool("hall.solid",false), true, true, true
	    );
    
    //----------------
    // the opening in the wall for the beam dump
    const double xCenterOnInWall = dump->enclosureCenterInMu2e()[0] 
      + tan(dump->coreRotY()) * (hallAirZmin - dump->enclosureCenterInMu2e()[2]);


    const double openingHalfWidth = dump->enclosureHalfSize()[0]/cos(dump->coreRotY());

    // the longer side of the "-x" slab
    const double snwdxLong = building->hallInsideHalfLengths()[0] + building->hallWallThickness()
      + (xCenterOnInWall - building->hallCenterInMu2e()[0]) - openingHalfWidth;

    // the shorter side of the "+x" slab
    const double spwdxShort = building->hallInsideHalfLengths()[0] + building->hallWallThickness()
      - (xCenterOnInWall - building->hallCenterInMu2e()[0]) - openingHalfWidth;

    // the difference between long and short sides
    const double sdx = building->hallWallThickness() * tan(dump->coreRotY());

    // "-x" slab
    VolumeInfo snx;
    snx.name = "HallNegZWallNXSlab";
    snx.solid = new G4Trap(snx.name,
			   2*dump->enclosureHalfSize()[1],
			   building->hallWallThickness(),
			   snwdxLong, snwdxLong - sdx);

    static CLHEP::HepRotation snxRotation(CLHEP::HepRotation::IDENTITY);
    snxRotation.rotateX(90*CLHEP::degree);

    finishNesting(snx, 
		  materialFinder.get("hall.wallMaterialName"),
		  &snxRotation,
		  CLHEP::Hep3Vector(0.5*((xCenterOnInWall - openingHalfWidth - 0.5*sdx) + 
					 (building->hallCenterInMu2e()[0]
					  - building->hallInsideHalfLengths()[0] - building->hallWallThickness())
					 ),
				    dump->enclosureCenterInMu2e()[1],
				    hallAirZmin - 0.5*building->hallWallThickness())
		  - parent.centerInMu2e(),
		  parent.logical,
		  0,
		  config.getBool("protonBeamDump.logicalEnclosureVisible"),
		  G4Colour::Red(),
		  false, true, true, true
		  );

    //----------------
    // "+x" slab
    VolumeInfo spx;
    spx.name = "HallNegZWallPXSlab";
    spx.solid = new G4Trap(spx.name,
			   2*dump->enclosureHalfSize()[1],
			   building->hallWallThickness(),
			   spwdxShort + sdx, spwdxShort);

    static CLHEP::HepRotation spxRotation(CLHEP::HepRotation::IDENTITY);
    spxRotation.rotateY(180*CLHEP::degree).rotateX(90*CLHEP::degree);

    finishNesting(spx, 
		  materialFinder.get("hall.wallMaterialName"),
		  &spxRotation,
		  CLHEP::Hep3Vector(0.5*((xCenterOnInWall + openingHalfWidth - 0.5*sdx) + 
					 (building->hallCenterInMu2e()[0]
					  + building->hallInsideHalfLengths()[0] + building->hallWallThickness())
					 ),
				    dump->enclosureCenterInMu2e()[1],
				    hallAirZmin - 0.5*building->hallWallThickness())
		  - parent.centerInMu2e(),
		  parent.logical,
		  0,
		  config.getBool("protonBeamDump.logicalEnclosureVisible"),
		  G4Colour::Red(),
		  false, true, true, true
		  );

    //----------------------------------------------------------------
    // Visual references for debugging
    const bool addTestObjects  = false;
    if(addTestObjects) {
      VolumeInfo ti;
      ti.name = "TestSphere1";
      ti.solid = new G4Orb(ti.name, 0.1*CLHEP::meter);
      finishNesting(ti, 
		    materialFinder.get("protonBeamDump.material.air"),
		    0,
		    CLHEP::Hep3Vector(xCenterOnInWall,
				      dump->enclosureCenterInMu2e()[1],
				      hallAirZmin - 0.5*building->hallWallThickness())
		    - parent.centerInMu2e(),
		    parent.logical, 0, true, G4Colour::Red(), true, true, true, true
		    );
    }
    if(addTestObjects) {
      VolumeInfo ti;
      ti.name = "TestSphere2";
      ti.solid = new G4Orb(ti.name, 0.1*CLHEP::meter);
      finishNesting(ti, 
		    materialFinder.get("protonBeamDump.material.air"),
		    0,
		    CLHEP::Hep3Vector(0, 0, 0),
		    wi.logical, 0, true, G4Colour::Blue(), true, true, true, true
		    );
    }
    if(addTestObjects) {
      VolumeInfo ti;
      ti.name = "TestCyl1";
      ti.solid = new G4Tubs(ti.name, 0., 5*CLHEP::cm, openingHalfWidth, 0., CLHEP::twopi);
      static CLHEP::HepRotation r(CLHEP::HepRotation::IDENTITY);
      r.rotateY(90*CLHEP::degree);
      
      finishNesting(ti, 
		    materialFinder.get("protonBeamDump.material.air"),
		    &r,
		    CLHEP::Hep3Vector(xCenterOnInWall,
				      dump->enclosureCenterInMu2e()[1],
				      hallAirZmin)
		    - parent.centerInMu2e(),
		    parent.logical, 0, true, G4Colour::Red(), true, true, true, true
		    );
    }
    if(addTestObjects) {
      VolumeInfo ti;
      ti.name = "TestCylPBDE";
      ti.solid = new G4Tubs(ti.name, 0., 5*CLHEP::cm, lcenter, 0., CLHEP::twopi);
      finishNesting(ti, 
		    materialFinder.get("protonBeamDump.material.air"),
		    0,
		    CLHEP::Hep3Vector(0, 0, 0),
		    logicalEnclosure.logical, 0, true, G4Colour::Blue(), true, true, true, true
		    );
    }
 
    //----------------------------------------------------------------
  }
}
