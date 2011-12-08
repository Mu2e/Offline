//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
// $Id: constructHall.cc,v 1.7 2011/12/08 17:20:58 gandr Exp $
// $Author: gandr $
// $Date: 2011/12/08 17:20:58 $
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
#include "GeometryService/inc/Mu2eBuilding.hh"
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

using namespace std;

namespace mu2e {

  VolumeInfo constructHall(const VolumeInfo& parent, SimpleConfig const * const config ) {
    
    bool const forceAuxEdgeVisible = config->getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = config->getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    MaterialFinder materialFinder(*config);

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
                                   config->getBool("hall.visible",true),
                                   G4Colour::Red(),
                                   config->getBool("hall.solid",false),
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   );



//tmp:    // Floor thickness.
//tmp:    const double ceilingThick = building->hallCeilingThickness();
//tmp:    const double floorThick   = building->hallFloorThickness();
//tmp:    const double wallThick    = building->hallWallThickness();
//tmp:
//tmp:    // Half lengths of the exterior of the concrete for the hall walls.
//tmp:    const double hallOutHLen[3] = {
//tmp:      hallFormalHLen[0] + wallThick,
//tmp:      hallFormalHLen[1] + (ceilingThick + floorThick)/2.,
//tmp:      hallFormalHLen[2] + wallThick /*the "-z" wall is not a part of hall concrete*/ - 0.5*wallThick
//tmp:    };
//tmp:
//tmp:    // Center of the concrete volume in the coordinate system of the parent.
//tmp:    const CLHEP::Hep3Vector concretePositionInParent = 
//tmp:      // Position of the center of the air in the world volume.
//tmp:      world->mu2eOriginInWorld() + building->hallCenterInMu2e()
//tmp:      // Correct for the possible asymmetry in the thickness of the floor/ceiling concrete
//tmp:      + CLHEP::Hep3Vector(0,  -(floorThick - ceilingThick)/2, 0)
//tmp:      // The wall on the negative Z side will be created by constructProtonBeamDump
//tmp:      // and is not a part of the hall concrete.
//tmp:      + CLHEP::Hep3Vector(0, 0, 0.5*wallThick)
//tmp:      // This makes it relative to the parent
//tmp:      - parent.centerInWorld
//tmp:      ;
//tmp:    
//tmp:    // Position of the hall air volume.  The only possible shift is 
//tmp:    // due to a non-equal thickness of the floor and the ceiling.
//tmp:    const CLHEP::Hep3Vector airPositionInConcrete(0, (floorThick - ceilingThick)/2, 
//tmp:						  // account for the missign "-z" wall
//tmp:						  -0.5*wallThick
//tmp:						  );
//tmp:
//tmp:    // Materials for the hall walls and the interior of the hall
//tmp:    G4Material* wallMaterial = materialFinder.get("hall.wallMaterialName");
//tmp:
//tmp:    // Concrete walls of the hall.
//tmp:    VolumeInfo wallInfo = nestBox( "HallWalls",
//tmp:                                   hallOutHLen,
//tmp:                                   wallMaterial,
//tmp:                                   0,
//tmp:                                   concretePositionInParent,
//tmp:                                   parent,
//tmp:                                   0,
//tmp:                                   hallVisible,
//tmp:                                   G4Colour::Red(),
//tmp:                                   hallSolid,
//tmp:                                   forceAuxEdgeVisible,
//tmp:                                   placePV,
//tmp:                                   doSurfaceCheck
//tmp:                                   );
//tmp:

//test:    //----------------------------------------------------------------
//test:    // TEST TEST TEST
//test:
//test:    static CLHEP::HepRotation testRotation(CLHEP::HepRotation::IDENTITY);
//test:    //testRotation.rotateX(-90*CLHEP::degree);
//test:
//test:    if(true) {
//test:      std::vector<G4TwoVector> s1;
//test:      s1.push_back(G4TwoVector(0, 0));
//test:      s1.push_back(G4TwoVector(0, 1000));
//test:      s1.push_back(G4TwoVector(1000, 0));
//test:      
//test:      VolumeInfo ti;
//test:      ti.name = "s1";
//test:      ti.solid = new G4ExtrudedSolid(ti.name, s1, 150, G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
//test:      
//test:      finishNesting(ti, 
//test:		    materialFinder.get("protonBeamDump.material.air"),
//test:		    &testRotation,
//test:		    CLHEP::Hep3Vector(3000, 0, 0),
//test:		    hallInfo.logical,
//test:		    0, true, G4Colour::Red(), true, true, true, true
//test:		    );
//test:      
//test:
//test:      const bool addTestObjects  = true;
//test:      if(addTestObjects) {
//test:	VolumeInfo si;
//test:	si.name = "TestSphere1";
//test:	si.solid = new G4Orb(si.name, 0.1*CLHEP::meter);
//test:	finishNesting(si, 
//test:		      materialFinder.get("protonBeamDump.material.air"),
//test:		      0,
//test:		      CLHEP::Hep3Vector(0, 0, 0),
//test:		      ti.logical, 0, true, G4Colour::Magenta(), true, true, true, true
//test:		      );
//test:      }
//test:    }
//test:    if(true) {
//test:      std::vector<G4TwoVector> s;
//test:      s.push_back(G4TwoVector(0, 0));
//test:      s.push_back(G4TwoVector(0, -1000));
//test:      s.push_back(G4TwoVector(-1000, 0));
//test:      
//test:      VolumeInfo ti;
//test:      ti.name = "s2";
//test:      ti.solid = new G4ExtrudedSolid(ti.name, s, 150, G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
//test:      
//test:      finishNesting(ti, 
//test:		    materialFinder.get("protonBeamDump.material.air"),
//test:		    &testRotation,
//test:		    CLHEP::Hep3Vector(3000, 0, 0),
//test:		    hallInfo.logical,
//test:		    0, true, G4Colour::Blue(), true, true, true, true
//test:		    );
//test:      
//test:    }
//test:    //----------------------------------------------------------------


    return hallInfo;

  }

}
