//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
// $Id: constructHall.cc,v 1.17 2013/09/27 13:23:26 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/09/27 13:23:26 $
//
// Original author KLG based on Mu2eWorld constructHall
//
// Notes:
// Construct the earthen overburden

// Mu2e includes
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "Mu2eG4/inc/constructHall.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Orb.hh"
#include "G4TwoVector.hh"
#include "CLHEP/Vector/Rotation.h"
#include "G4NistManager.hh"

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
    G4Material* ceilingMaterial = materialFinder.get("hall.ceilingMaterialName");

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

    std::vector<G4TwoVector> horizontalConcreteOutlineExt; // ceiling extension
    std::copy(building->concreteOuterOutline1().begin(),
              building->concreteOuterOutline1().end(),
              std::back_inserter(horizontalConcreteOutlineExt));

    std::copy(building->concreteOuterOutlineExt().begin(),
              building->concreteOuterOutlineExt().end(),
              std::back_inserter(horizontalConcreteOutlineExt));

    std::copy(building->concreteOuterOutline3().begin(),
              building->concreteOuterOutline3().end(),
              std::back_inserter(horizontalConcreteOutlineExt));

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
                  ceilingMaterial,
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

    VolumeInfo hallCeilingExt("HallConcreteCeilingExt",
                           CLHEP::Hep3Vector(0, building->hallInsideYmax() + 3*building->hallCeilingThickness()/2, 0)
                           - hallInfo.centerInMu2e(),
                           hallInfo.centerInWorld);

    hallCeilingExt.solid = new G4ExtrudedSolid(hallCeilingExt.name,
                                               horizontalConcreteOutlineExt,
                                               building->hallCeilingThickness()/2,
                                               G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

    finishNesting(hallCeilingExt,
                  ceilingMaterial,
                  &horizontalConcreteRotation,
                  hallCeilingExt.centerInParent,
                  hallInfo.logical,
                  0,
                  config.getBool("hall.ceilingVisible"),
                  G4Colour::Grey(),
                  config.getBool("hall.ceilingSolid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );


    // Create beamline shielding slabs
    //
    // - EMFB dirt overlaps first slab.  Introduce kludge to make
    // - first slab not overlab with dirt; and the remaining part of
    // - the first slab extend the full intended width (36')
    double yElev = building->hallInsideYmax() + 2*building->hallCeilingThickness();

    // Get height of dirt
    GeomHandle<ProtonBeamDump> dump;
    GeomHandle<ExtMonFNALBuilding> emfb;
    //    GeomHandle<WorldG4> world;
    
    //    const double dumpDirtYmin =
    //      world->hallFormalCenterInWorld()[1] - world->hallFormalHalfSize()[1]
    //      - world->mu2eOriginInWorld()[1]
      ;
    const double dumpDirtYmax = std::max(
                                         emfb->roomInsideYmax() + emfb->roomCeilingThickness() + emfb->dirtOverheadThickness()
                                         ,
                                         dump->frontShieldingCenterInMu2e()[1] + dump->frontShieldingHalfSize()[1]
                                         );
    
    const Box &  slab2                   = building->concreteBeamlineSlab( 1 );
    const double slab2TopHeight          = yElev + 4*slab2.getYhalfLength();
    const double dumpDirtSlab2HeightDiff = slab2TopHeight-dumpDirtYmax;

    // Portion up to dirt level

    // Portion above dirt level

    Box box( slab2.getXhalfLength(), 0.5*dumpDirtSlab2HeightDiff, slab2.getZhalfLength() );
    CLHEP::Hep3Vector pos = CLHEP::Hep3Vector(building->xPosOfSlabEnd()-box.getXhalfLength(), dumpDirtYmax+0.5*dumpDirtSlab2HeightDiff, 0)-hallInfo.centerInMu2e();
    
    cout << " Slab 2b: " << pos << "   " << box.getXhalfLength() << " " << box.getYhalfLength() << " " << box.getZhalfLength() << endl;

    nestBox( "BeamlineSlab_2b",
             box,
             ceilingMaterial,
             0,
             pos, // position wrt hall
             hallInfo,
             0,
             config.getBool("hall.ceilingVisible"),
             G4Colour::Grey(),
             config.getBool("hall.ceilingSolid"),
             forceAuxEdgeVisible,
             placePV,
             doSurfaceCheck );
    
    yElev += 4*slab2.getYhalfLength();


    for ( std::size_t iSlab = 2 ; iSlab < building->nBeamlineSlabs() ; iSlab++ ) {
      ostringstream slab;
      slab << "BeamlineSlab" << iSlab+1;

      Box const& box = building->concreteBeamlineSlab( iSlab );

      yElev += box.getYhalfLength();
      CLHEP::Hep3Vector pos = CLHEP::Hep3Vector(building->xPosOfSlabEnd()-box.getXhalfLength(), yElev, 0)-hallInfo.centerInMu2e();

      cout << slab.str() << ": " << pos << "   " << box.getXhalfLength() << " " << box.getYhalfLength() << " " << box.getZhalfLength() << endl;

      nestBox( slab.str(),
               box,
               ceilingMaterial,
               0,
               pos, // position wrt hall
               hallInfo,
               0,
               config.getBool("hall.ceilingVisible"),
               G4Colour::Grey(),
               config.getBool("hall.ceilingSolid"),
               forceAuxEdgeVisible,
               placePV,
               doSurfaceCheck );

      yElev += box.getYhalfLength();
    }

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
