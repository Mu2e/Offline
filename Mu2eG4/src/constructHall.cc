//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
// $Id: constructHall.cc,v 1.18 2013/09/27 17:19:39 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/09/27 17:19:39 $
//
// Original author KLG based on Mu2eWorld constructHall
//
// Notes:
// Construct the earthen overburden

// Mu2e includes
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eBuildingGeom/inc/Mu2eHall.hh"
#include "Mu2eG4/inc/constructHall.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/nestBox.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Orb.hh"
#include "G4TwoVector.hh"
#include "CLHEP/Vector/Rotation.h"
#include "G4NistManager.hh"

// C++ includes
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

namespace mu2e {

  VolumeInfo constructHall(const VolumeInfo& worldInfo, const SimpleConfig& config ) {

    MaterialFinder materialFinder( config );

    GeomHandle<WorldG4> world;
    GeomHandle<Mu2eHall> building;

    const auto& geoOptions = art::ServiceHandle<GeometryService>()->geomOptions();

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
                                   geoOptions->forceAuxEdgeVisible( "HallAir" ),
                                   geoOptions->placePV( "HallAir" ),
                                   geoOptions->doSurfaceCheck( "HallAir" )
                                   );

    // Rotation is static because rotations are not copied into G4.
    static CLHEP::HepRotation horizontalConcreteRotation(CLHEP::HepRotation::IDENTITY);
    horizontalConcreteRotation.rotateX( 90*CLHEP::degree);
    horizontalConcreteRotation.rotateZ( 90*CLHEP::degree);

    constructSolids( hallInfo, building->getBldgSolids(), horizontalConcreteRotation );
    constructSolids( hallInfo, building->getDirtSolids(), horizontalConcreteRotation );

    return hallInfo;

  }

  //================================================================================
  void constructSolids( const VolumeInfo& hallInfo, 
			const std::map<std::string,ExtrudedSolid>& solidMap,
			const CLHEP::HepRotation& rot) {
    
    //-----------------------------------------------------------------
    // Building and dirt volumes are extruded solids.
    //-----------------------------------------------------------------
    
    const auto& geoOptions = art::ServiceHandle<GeometryService>()->geomOptions();

    for ( const auto& keyVolumePair : solidMap ) {

      const auto& volume = keyVolumePair.second;

      VolumeInfo tmpVol(volume.getName(),
                        volume.getOffsetFromMu2eOrigin() - hallInfo.centerInMu2e(),
                        hallInfo.centerInWorld);
      
      tmpVol.solid = new G4ExtrudedSolid(tmpVol.name, 
                                         volume.getVertices(),
                                         volume.getYhalfThickness(),
                                         G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
      
      finishNesting(tmpVol,
                    findMaterialOrThrow( volume.getMaterial() ),
                    &rot,
                    tmpVol.centerInParent,
                    hallInfo.logical,
                    0,
                    geoOptions->isVisible( volume.getName() ),
                    G4Colour::Grey(),
                    geoOptions->isSolid( volume.getName() ),
                    geoOptions->forceAuxEdgeVisible( volume.getName() ),
                    geoOptions->placePV( volume.getName() ),
                    geoOptions->doSurfaceCheck( volume.getName() )
                    );
    }

  }

}
