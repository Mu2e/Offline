// Constructs dirt overburden inside the formal hall box.
// Note that there is also dirt around the hall box.
//
// $Id: constructDirt.cc,v 1.10 2013/09/27 20:59:34 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/09/27 20:59:34 $
//
// Original author KLG based on Mu2eWorld constructDirt
// Updated by Andrei Gaponenko.

// Mu2e includes.
#include "Mu2eG4/inc/constructDirt.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Paraboloid.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;

namespace mu2e {

  void constructDirt(const VolumeInfo& parent, const SimpleConfig& config) {

    G4Helper* _helper = &(*art::ServiceHandle<G4Helper>() );

    // Here we can e.g. place dirt on top of the hall ceiling.
    GeomHandle<WorldG4> world;
    GeomHandle<Mu2eBuilding> building;
    GeomHandle<ProtonBeamDump> dump;
    GeomHandle<ExtMonFNALBuilding> emfb;

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;
    
    const double xmax = world->halfLengths()[0];
    const double xmin = -world->hallFormalCenterInWorld().x() + world->hallFormalHalfSize()[0];
    const double deltaX = xmax-xmin;
    const double zmin = -world->halfLengths()[2];
    const double zmax = world->hallFormalCenterInWorld().z() - world->hallFormalHalfSize()[2];
    const double deltaZ = zmax-zmin;

    const double zExtentBuilding = building->hallInsideZDSCorner() + building->hallWallThickness();
    const double yExtentDirt     = building->hallInsideYmax() + 2*building->hallCeilingThickness();
    const double xExtentBuilding = building->hallInsideXDSCorner() + building->hallWallThickness();
    const double zExtent    =  building->hallInsideZmax()+building->hallWallThickness();
    const double yExtentLow =  building->hallInsideYmin();
    const double xExtent    =  building->hallInsideXmax()+building->hallWallThickness();
    
    

    const double dirtDimensions[3] = { 0.5*(xExtent-xExtentBuilding),
                                       0.5*(yExtentDirt - yExtentLow),
                                       0.5*(zExtent - zExtentBuilding) };

    CLHEP::Hep3Vector loc( xExtentBuilding+dirtDimensions[0], 
                           yExtentLow+dirtDimensions[1], 
                           zExtentBuilding+dirtDimensions[2] );
    
    MaterialFinder materialFinder(config);

    // The formal hall volume
    VolumeInfo hallInfo = nestBox( "CRV_dirt",
                                   dirtDimensions,
                                   materialFinder.get("dirt.overburdenMaterialName"),
                                   0,
                                   loc-_helper->locateVolInfo("HallAir").centerInMu2e(),
                                   _helper->locateVolInfo("HallAir"),
                                   0,
                                   true,
                                   G4Colour::Brown(),
                                   false,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   );

  } // constructDirt()

} // namespace mu2e
