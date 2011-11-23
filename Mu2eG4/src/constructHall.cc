//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
// $Id: constructHall.cc,v 1.6 2011/11/23 16:42:20 gandr Exp $
// $Author: gandr $
// $Date: 2011/11/23 16:42:20 $
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

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"

using namespace std;

namespace mu2e {

  VolumeInfo constructHall( const VolumeInfo& dirt,
                            SimpleConfig const * const _config
                            ){

    // A helper class.
    MaterialFinder materialFinder(*_config);

    GeomHandle<WorldG4> world;
    GeomHandle<Mu2eBuilding> building;

    const vector<double>& hallInHLen = building->hallInsideHalfLengths();

    // Floor thickness.
    const double ceilingThick = building->hallCeilingThickness();
    const double floorThick   = building->hallFloorThickness();
    const double wallThick    = building->hallWallThickness();

    // Half lengths of the exterior of the concrete for the hall walls.
    const double hallOutHLen[3] = {
      hallInHLen[0] + wallThick,
      hallInHLen[1] + (ceilingThick + floorThick)/2.,
      hallInHLen[2] + wallThick /*the "-z" wall is not a part of hall concrete*/ - 0.5*wallThick
    };

    // Center of the concrete volume in the coordinate system of the dirt.
    const CLHEP::Hep3Vector concretePositionInDirt = 
      // Position of the center of the air in the world volume.
      world->mu2eOriginInWorld() + building->hallCenterInMu2e()
      // Correct for the possible asymmetry in the thickness of the floor/ceiling concrete
      + CLHEP::Hep3Vector(0,  -(floorThick - ceilingThick)/2, 0)
      // The wall on the negative Z side will be created by constructProtonBeamDump
      // and is not a part of the hall concrete.
      + CLHEP::Hep3Vector(0, 0, 0.5*wallThick)
      // This makes it relative to the dirt
      - dirt.centerInWorld
      ;
    
    // Position of the hall air volume.  The only possible shift is 
    // due to a non-equal thickness of the floor and the ceiling.
    const CLHEP::Hep3Vector airPositionInConcrete(0, (floorThick - ceilingThick)/2, 
						  // account for the missign "-z" wall
						  -0.5*wallThick
						  );

    // Materials for the hall walls and the interior of the hall
    G4Material* wallMaterial = materialFinder.get("hall.wallMaterialName");
    G4Material* hallMaterial = materialFinder.get("hall.insideMaterialName");

    bool hallVisible = _config->getBool("hall.visible",true);
    bool hallSolid   = _config->getBool("hall.solid",false);

    bool const forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    // Concrete walls of the hall.
    VolumeInfo wallInfo = nestBox( "HallWalls",
                                   hallOutHLen,
                                   wallMaterial,
                                   0,
                                   concretePositionInDirt,
                                   dirt,
                                   0,
                                   hallVisible,
                                   G4Colour::Red(),
                                   hallSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   );

    // Air volume inside of the hall.
    VolumeInfo hallInfo = nestBox( "HallAir",
                                   hallInHLen,
                                   hallMaterial,
                                   0,
                                   airPositionInConcrete,
                                   wallInfo,
                                   0,
                                   hallVisible,
                                   G4Colour::Red(),
                                   hallSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   );

    return hallInfo;

  }

}
