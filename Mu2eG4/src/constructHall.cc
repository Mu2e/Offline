//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
// $Id: constructHall.cc,v 1.13 2012/04/17 19:56:56 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/17 19:56:56 $
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

  VolumeInfo constructHall(const VolumeInfo& worldInfo, const SimpleConfig& config) {

    bool const forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    bool const doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    bool const placePV             = true;

    MaterialFinder materialFinder(config);

    // Materials for the hall walls, floor, and ceiling
    G4Material* wallMaterial = materialFinder.get("hall.wallMaterialName");

    GeomHandle<WorldG4> world;
    GeomHandle<Mu2eBuilding> building;

    // The formal hall volume
    VolumeInfo hallInfo = nestBox( "HallAir",
                                   world->hallFormalHalfSize(),
                                   materialFinder.get("hall.insideMaterialName"),
                                   0,
                                   world->hallFormalCenterInWorld(),
                                   worldInfo,
                                   0,
                                   config.getBool("hall.formalBoxVisible"),
                                   G4Colour::Red(),
                                   config.getBool("hall.formalBoxSolid"),
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
                         - hallInfo.centerInMu2e(),
                         hallInfo.centerInWorld);

    hallFloor.solid = new G4ExtrudedSolid(hallFloor.name, horizontalConcreteOutline,
                                          building->hallFloorThickness()/2,
                                          G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

    finishNesting(hallFloor,
                  wallMaterial,
                  &horizontalConcreteRotation,
                  hallFloor.centerInParent,
                  hallInfo.logical,
                  0,
                  config.getBool("hall.floorVisible"),
                  G4Colour::Grey(),
                  config.getBool("hall.floorSolid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );


    VolumeInfo hallCeiling("HallConcreteCeiling",
                           CLHEP::Hep3Vector(0, building->hallInsideYmax() + building->hallCeilingThickness()/2, 0)
                           - hallInfo.centerInMu2e(),
                           hallInfo.centerInWorld);

    hallCeiling.solid = new G4ExtrudedSolid(hallCeiling.name, horizontalConcreteOutline,
                                            building->hallCeilingThickness()/2,
                                            G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

    finishNesting(hallCeiling,
                  wallMaterial,
                  &horizontalConcreteRotation,
                  hallCeiling.centerInParent,
                  hallInfo.logical,
                  0,
                  config.getBool("hall.ceilingVisible"),
                  G4Colour::Grey(),
                  config.getBool("hall.ceilingSolid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    //----------------------------------------------------------------
    // Create hall walls. The walls are placed inside the formal hall box.

    std::vector<G4TwoVector> wallConcreteOutline(building->hallInsideOutline());

    // Need the points go clockwise w.r.t. the wall concrete.
    std::reverse(wallConcreteOutline.begin(), wallConcreteOutline.end());

    // Remove the duplicate
    if(wallConcreteOutline.back() == horizontalConcreteOutline.front()) {
      wallConcreteOutline.pop_back();
    }

    std::copy(horizontalConcreteOutline.begin(),
              horizontalConcreteOutline.end(),
              std::back_inserter(wallConcreteOutline));

    static CLHEP::HepRotation hallWallsRotation(CLHEP::HepRotation::IDENTITY);
    hallWallsRotation.rotateX(-90*CLHEP::degree);

    VolumeInfo hallWalls("HallWalls",
                         CLHEP::Hep3Vector(0, (building->hallInsideYmin() + building->hallInsideYmax())/2, 0)
                         - hallInfo.centerInMu2e(),
                            hallInfo.centerInWorld);

    hallWalls.solid = new G4ExtrudedSolid(hallWalls.name, wallConcreteOutline,
                                             (building->hallInsideYmax() - building->hallInsideYmin())/2,
                                             G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

    finishNesting(hallWalls,
                  materialFinder.get("hall.wallMaterialName"),
                  &hallWallsRotation,
                  hallWalls.centerInParent,
                  hallInfo.logical,
                  0,
                  config.getBool("hall.wallsVisible"),
                  G4Colour(0.8, 0.8, 0.8), // lighter grey
                  config.getBool("hall.wallsSolid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    //----------------------------------------------------------------
    return hallInfo;
  }

}
