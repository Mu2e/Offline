//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
// $Id: constructHall.cc,v 1.2 2011/05/17 15:36:01 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:01 $
//
// Original author KLG based on Mu2eWorld constructHall
//
// Notes:
// Construct the earthen overburden

// Mu2e includes.
#include "Mu2eG4/inc/constructHall.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestBox.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"

using namespace std;

namespace mu2e {

  VolumeInfo constructHall( const VolumeInfo& parent, 
                            SimpleConfig const * const _config
                            ){

    // A helper class.
    MaterialFinder materialFinder(*_config);

    // Dimensions of the world.
    vector<double> worldHLen;
    _config->getVectorDouble("world.halfLengths", worldHLen, 3);

    vector<double> hallInHLen;
    _config->getVectorDouble("hall.insideHalfLengths",hallInHLen,3);

    // Floor thickness.
    double ceilingThick = _config->getDouble("hall.ceilingThick");
    double floorThick   = _config->getDouble("hall.floorThick");
    double wallThick    = _config->getDouble("hall.wallThick");

    // Top of the floor in G4 world coordinates.
    double yFloor = -worldHLen[1] + floorThick;

    // Position of the center of the hall in the world volume.
    vector<double> hallPosition;
    _config->getVectorDouble("hall.offset", hallPosition,3);
    double hallY0 = yFloor + hallInHLen[1] + hallPosition[1];

    // Materials for the hall walls and the interior of the hall
    G4Material* wallMaterial = materialFinder.get("hall.wallMaterialName");
    G4Material* hallMaterial = materialFinder.get("hall.insideMaterialName");

    // Half lengths of the exterior of the concrete for the hall walls.
    double hallOutHLen[3] ={
      hallInHLen[0] + wallThick,
      hallInHLen[1] + ( ceilingThick + floorThick )/2.,
      hallInHLen[2] + wallThick
    };
    
    // Center of the concrete volume in the coordinate system of the dirt.
    G4ThreeVector wallOffset = 
      G4ThreeVector(hallPosition[0], hallY0, hallPosition[2]) - parent.centerInParent;

    // Origin of the hall air volume in the system of the hall concrete volume.
    G4ThreeVector hallOffset( 0., (floorThick-ceilingThick)/2., 0.);

    bool hallVisible = _config->get<bool>("hall.visible",true);
    bool hallSolid   = _config->get<bool>("hall.solid",false);

    bool const forceAuxEdgeVisible = _config->get<bool>("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config->get<bool>("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    // Concrete walls of the hall.
    VolumeInfo wallInfo = nestBox( "HallWalls",
                                   hallOutHLen,
                                   wallMaterial,
                                   0,
                                   wallOffset,
                                   parent,
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
                                   hallOffset,
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
