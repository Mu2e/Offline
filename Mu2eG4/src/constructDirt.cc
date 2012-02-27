//
// Free function to create the earthen overburden.
//
// $Id: constructDirt.cc,v 1.8 2012/02/27 06:05:35 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/27 06:05:35 $
//
// Original author KLG based on Mu2eWorld constructDirt
//
// Notes:
//  The Earth overburden is modeled in two parts: a box that extends
//  to the surface of the earth plus a cap above grade.  The cap is shaped
//  as a G4Paraboloid.
//
// Mu2e includes.
#include "Mu2eG4/inc/constructDirt.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Paraboloid.hh"

using namespace std;

namespace mu2e {

  VolumeInfo constructDirt( const VolumeInfo& parent,
                            SimpleConfig const * const _config
                            ){

    // Dimensions and material of the world.
    GeomHandle<WorldG4> world;
    GeomHandle<Mu2eBuilding> building;
    G4Material* dirtMaterial = MaterialFinder(*_config).get("dirt.overburdenMaterialName");

    const bool dirtVisible    = _config->getBool("dirt.visible",true);
    const bool dirtSolid      = _config->getBool("dirt.solid",false);
    const bool dirtCapVisible = _config->getBool("dirt.capVisible",true);
    const bool dirtCapSolid   = _config->getBool("dirt.capSolid",false);
    const bool forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    const bool placePV             = true;

    // dirt box depth, to the bottom of the world
    const double dirtTotalDepth = world->halfLengths()[1] + world->dirtG4Ymax();

    // Half lengths of the dirt box.
    double dirtHLen[3] = { world->halfLengths()[0], 0.5*dirtTotalDepth, world->halfLengths()[2] };

    // Center of the dirt box, in the G4 world system.
    G4ThreeVector dirtOffset(0., -world->halfLengths()[1] + 0.5*dirtTotalDepth, 0.);

    // Main body of dirt around the hall.
    VolumeInfo dirtInfo = nestBox( "DirtBody",
                                   dirtHLen,
                                   dirtMaterial,
                                   0,
                                   dirtOffset,
                                   parent,
                                   0,
                                   dirtVisible,
                                   G4Colour::Magenta(),
                                   dirtSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   );

    // Dirt cap is modeled as a paraboloid.
    double capHalfHeight = building->dirtCapHalfHeight();
    double capBottomR    = building->dirtCapBottomRadius();
    double capTopR       = building->dirtCapTopRadius();

    double dsz0          = _config->getDouble("toyDS.z0");

    GeomHandle<Beamline> beamg;
    double solenoidOffset = beamg->solenoidOffset();

    G4ThreeVector dirtCapOffset( -solenoidOffset,
				 -world->halfLengths()[1] + dirtTotalDepth + capHalfHeight,
                                 dsz0 + world->mu2eOriginInWorld().z());

    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
    G4RotationMatrix* dirtCapRot = reg.add(G4RotationMatrix());
    dirtCapRot->rotateX( -90*CLHEP::degree);

    string dirtCapName("DirtCap");

    // Construct the cap.
    VolumeInfo dirtCapInfo( dirtCapName, dirtCapOffset, dirtInfo.centerInWorld);

    dirtCapInfo.solid = new G4Paraboloid( dirtCapName, capHalfHeight, capTopR, capBottomR);

    finishNesting(dirtCapInfo,
                  dirtMaterial,
                  dirtCapRot,
                  dirtCapOffset,
                  parent.logical,
                  0,
                  dirtCapVisible,
                  G4Colour::Green(),
                  dirtCapSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    return dirtInfo;

  }

}
