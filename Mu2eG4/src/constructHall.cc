//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
// $Id: constructHall.cc,v 1.9 2012/02/27 06:05:35 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/27 06:05:35 $
//
// Original author KLG based on Mu2eWorld constructHall
//
// Notes:
// Construct the earthen overburden

// Mu2e includes.
#include "Mu2eG4/inc/constructHall.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Orb.hh"
#include "G4TwoVector.hh"
#include "CLHEP/Vector/Rotation.h"

#include <vector>
#include <algorithm>

using namespace std;

namespace mu2e {

  VolumeInfo constructHall(const VolumeInfo& parent, SimpleConfig const * const config ) {
    
    bool const forceAuxEdgeVisible = config->getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = config->getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    MaterialFinder materialFinder(*config);

    // Materials for the hall walls, floor, and ceiling
    G4Material* wallMaterial = materialFinder.get("hall.wallMaterialName");

    GeomHandle<WorldG4> world;
    GeomHandle<Mu2eBuilding> building;

    vector<double> hallFormalHLen(3);
    hallFormalHLen[0] = (building->hallInsideXmax() - building->hallInsideXmin())/2;
    hallFormalHLen[1] = (building->hallInsideYmax() - building->hallInsideYmin())/2;
    hallFormalHLen[2] = (building->hallInsideZmax() - world->hallFormalZminInMu2e())/2;

    CLHEP::Hep3Vector hallFormalCenter
      ( (building->hallInsideXmax() + building->hallInsideXmin())/2,
	(building->hallInsideYmax() + building->hallInsideYmin())/2,
	(building->hallInsideZmax() + world->hallFormalZminInMu2e())/2
	);

    // The formal hall volume
    VolumeInfo hallInfo = nestBox( "HallAir",
                                   hallFormalHLen,
                                   materialFinder.get("hall.insideMaterialName"),
                                   0,
                                   hallFormalCenter - parent.centerInMu2e(),
                                   parent,
                                   0,
                                   config->getBool("hall.formalBoxVisible",false),
                                   G4Colour::Red(),
                                   config->getBool("hall.formalBoxSolid",false),
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   );

    //----------------------------------------------------------------
    // Floor and ceiling concrete are extruded solids.
    // The concrete extends beyond the hall air outline to "support" the walls.

    // points need to be in the clock-wise order
    std::vector<G4TwoVector> horizontalConcreteOutline; // same for floor and ceiling
    std::copy(building->concreteOuterOutline1().begin(), 
	      building->concreteOuterOutline1().end(),
	      std::back_inserter(horizontalConcreteOutline));

    std::copy(building->concreteOuterOutline2().begin(), 
	      building->concreteOuterOutline2().end(),
	      std::back_inserter(horizontalConcreteOutline));

    std::copy(building->concreteOuterOutline3().begin(), 
	      building->concreteOuterOutline3().end(),
	      std::back_inserter(horizontalConcreteOutline));

    //----------------
    static CLHEP::HepRotation horizontalConcreteRotation(CLHEP::HepRotation::IDENTITY);
    horizontalConcreteRotation.rotateX(-90*CLHEP::degree);
    
    VolumeInfo hallFloor("HallConcreteFloor",
			 CLHEP::Hep3Vector(0, building->hallInsideYmin() - building->hallFloorThickness()/2, 0)
			 - parent.centerInMu2e(),
			 parent.centerInWorld);

    hallFloor.solid = new G4ExtrudedSolid(hallFloor.name, horizontalConcreteOutline, 
					  building->hallFloorThickness()/2, 
					  G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
    
    finishNesting(hallFloor,
		  wallMaterial,
		  &horizontalConcreteRotation,
		  hallFloor.centerInParent,
		  parent.logical,
		  0, 
		  config->getBool("hall.floorVisible",true),
		  G4Colour::Grey(),
		  config->getBool("hall.floorSolid",false),
		  forceAuxEdgeVisible,
		  placePV,
		  doSurfaceCheck
		  );


    VolumeInfo hallCeiling("HallConcreteCeiling",
			   CLHEP::Hep3Vector(0, building->hallInsideYmax() + building->hallCeilingThickness()/2, 0)
			   - parent.centerInMu2e(),
			   parent.centerInWorld);
    
    hallCeiling.solid = new G4ExtrudedSolid(hallCeiling.name, horizontalConcreteOutline, 
					    building->hallCeilingThickness()/2, 
					    G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
    
    finishNesting(hallCeiling,
		  wallMaterial,
		  &horizontalConcreteRotation,
		  hallCeiling.centerInParent,
		  parent.logical,
		  0, 
		  config->getBool("hall.ceilingVisible",true),
		  G4Colour::Grey(),
		  config->getBool("hall.ceilingSolid",false),
		  forceAuxEdgeVisible,
		  placePV,
		  doSurfaceCheck
		  );

    //----------------------------------------------------------------
    // The two walls parallel to the YZ plane are simple boxes that extend
    // in Z outside of hall air to overlap with the XY walls.

    if(true) {
      std::vector<double> hl(3);
      hl[0] = building->hallWallThickness()/2;
      hl[1] = (building->hallInsideYmax() - building->hallInsideYmin())/2;
      hl[2] = (building->hallInsideZmax() - building->hallInsideZExtMonUCIWall())/2 + building->hallWallThickness();

      nestBox("HallConcreteWallPlusX",
	      hl,
	      wallMaterial,
	      0,
	      CLHEP::Hep3Vector(building->hallInsideXmax() + building->hallWallThickness()/2,
				(building->hallInsideYmax() + building->hallInsideYmin())/2,
				(building->hallInsideZmax() + building->hallInsideZExtMonUCIWall())/2
				)
	      - parent.centerInMu2e(),
	      parent,
	      0,
	      config->getBool("hall.wallsVisible", true),
	      G4Colour::Grey(),
	      config->getBool("hall.wallsSolid", false),
	      forceAuxEdgeVisible,
	      placePV,
	      doSurfaceCheck
	      );
    }


    if(true) {
      std::vector<double> hl(3);
      hl[0] = building->hallWallThickness()/2;
      hl[1] = (building->hallInsideYmax() - building->hallInsideYmin())/2;
      hl[2] = (building->hallInsideZmax() - building->hallInsideZBeamDumpWall())/2 + building->hallWallThickness();

      nestBox("HallConcreteWallMinusX",
	      hl,
	      wallMaterial,
	      0,
	      CLHEP::Hep3Vector(building->hallInsideXmin() - building->hallWallThickness()/2,
				(building->hallInsideYmax() + building->hallInsideYmin())/2,
				(building->hallInsideZmax() + building->hallInsideZBeamDumpWall())/2
				)
	      - parent.centerInMu2e(),
	      parent,
	      0,
	      config->getBool("hall.wallsVisible", true),
	      G4Colour::Grey(),
	      config->getBool("hall.wallsSolid", false),
	      forceAuxEdgeVisible,
	      placePV,
	      doSurfaceCheck
	      );
    }

    //----------------------------------------------------------------
    // The wall at Zmax is a simple box that fits in between floor,
    // ceiling, and and other walls.
    // Walls on the opposite side of hall are handled by the beam dump code.

    if(true) {
      std::vector<double> hl(3);
      hl[0] = (building->hallInsideXmax() - building->hallInsideXmin())/2;
      hl[1] = (building->hallInsideYmax() - building->hallInsideYmin())/2;
      hl[2] = building->hallWallThickness()/2;

      nestBox("HallConcreteWallPlusZ",
	      hl,
	      wallMaterial,
	      0,
	      CLHEP::Hep3Vector((building->hallInsideXmax() + building->hallInsideXmin())/2,
				(building->hallInsideYmax() + building->hallInsideYmin())/2,
				building->hallInsideZmax() + building->hallWallThickness()/2
				)
	      - parent.centerInMu2e(),
	      parent,
	      0,
	      config->getBool("hall.wallsVisible", true),
	      G4Colour::Grey(),
	      config->getBool("hall.wallsSolid", false),
	      forceAuxEdgeVisible,
	      placePV,
	      doSurfaceCheck
	      );
    }

    //----------------------------------------------------------------
    return hallInfo;
  }

}
