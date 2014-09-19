// Constructs CRV dirt burden inside the formal hall box, filling in
// the space +z wrt the TS left region, and +x wrt the DS left
// region. Note that there is also dirt around the hall box.
//
// $Id: constructDirt.cc,v 1.12 2014/09/19 19:15:00 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/19 19:15:00 $
//
// Original author KLG based on Mu2eWorld constructDirt
// Updated by Kyle Knoepfel

// Mu2e includes.
#include "Mu2eG4/inc/constructDirt.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"

// G4 includes
#include "G4Color.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;

namespace mu2e {

  void constructDirt(const VolumeInfo& parent, const SimpleConfig& config) {

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");

    const bool visible             = config.getBool("world.dirt.isVisible");
    const bool solid               = config.getBool("world.dirt.isSolid");
    const bool placePV             = true;

    // Get extent of CRV dirt
    GeomHandle<Mu2eBuilding> building;

    const double xExtentBuilding =  building->hallInsideXDSCorner() + building->hallWallThickness();
    const double xExtent         =  building->hallInsideXmax()+building->hallWallThickness();
    
    const double yExtentDirt     =  building->hallInsideYmax() + 2*building->hallCeilingThickness();
    const double yExtentLow      =  building->hallInsideYmin() - building->hallFloorThickness();

    const double zExtentBuilding =  building->hallInsideZDSCorner() + building->hallWallThickness();
    const double zExtent         =  building->hallInsideZmax()+building->hallWallThickness();
    
    const double dirtDimensions[3] = { 0.5*(xExtent-xExtentBuilding),
                                       0.5*(yExtentDirt - yExtentLow),
                                       0.5*(zExtent - zExtentBuilding) };

    const CLHEP::Hep3Vector loc( xExtentBuilding+dirtDimensions[0], 
                                 yExtentLow     +dirtDimensions[1], 
                                 zExtentBuilding+dirtDimensions[2] );
    
    art::ServiceHandle<GeometryService>()->geomOptions()->loadEntry( config, "CRV_dirt", "world.dirt" );

    nestBox( "CRV_dirt",
             dirtDimensions,
             findMaterialOrThrow( config.getString("dirt.overburdenMaterialName") ),
             0,
             loc-parent.centerInMu2e(),
             parent,
             0,
             visible,
             G4Colour::Brown(),
             solid,
             forceAuxEdgeVisible,
             placePV,
             doSurfaceCheck
             );

  } // constructDirt()

} // namespace mu2e
