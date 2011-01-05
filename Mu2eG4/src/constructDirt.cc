//
// Free function to create the earthen overburden.
//
// $Id: constructDirt.cc,v 1.1 2011/01/05 21:04:47 genser Exp $
// $Author: genser $
// $Date: 2011/01/05 21:04:47 $
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
#include "GeometryService/inc/GeometryService.hh"
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
    vector<double> worldHLen;
    _config->getVectorDouble("world.halfLengths", worldHLen, 3);

    // A helper class.
    MaterialFinder materialFinder(*_config);

    bool const forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    // Get parameters related to the overall dimensions of the hall and to
    // the earthen overburden.
    double floorThick           = CLHEP::mm * _config->getDouble("hall.floorThick");
    double ceilingThick         = CLHEP::mm * _config->getDouble("hall.ceilingThick");
    //double wallThick            = CLHEP::mm * _config->getDouble("hall.wallThick");
    double overburdenDepth      = CLHEP::mm * _config->getDouble("dirt.overburdenDepth");
    vector<double> hallInHLen;
    _config->getVectorDouble("hall.insideHalfLengths",hallInHLen,3);

    // Derived parameters.
    G4Material* dirtMaterial = materialFinder.get("dirt.overburdenMaterialName");

    // Top of the floor in G4 world coordinates.
    double yFloor   = -worldHLen[1] + floorThick;

    // The height above the floor of the y origin of the Mu2e coordinate system.
    //    double yOriginHeight = _config->getDouble("world.mu2eOrigin.height" )*CLHEP::mm;

    // Bottom of the ceiling in G4 world coordinates.
    double yCeilingInSide = yFloor + 2.*hallInHLen[1];
    
    // Top of the ceiling in G4 world coordinates.
    double yCeilingOutside  = yCeilingInSide + ceilingThick;

    // Surface of the earth in G4 world coordinates.
    double ySurface  = yCeilingOutside + overburdenDepth;
    
    // Half length and y origin of the dirt box.
    double yLDirt = ( ySurface + worldHLen[1] )/2.;
    double y0Dirt = -worldHLen[1] + yLDirt;
    
    // Center of the dirt box, in the G4 world system.
    G4ThreeVector dirtOffset(0.,y0Dirt,0.);
    
    // Half lengths of the dirt box.
    double dirtHLen[3] = { worldHLen[0], yLDirt, worldHLen[2] };

    bool dirtVisible    = _config->getBool("dirt.visible",true);
    bool dirtSolid      = _config->getBool("dirt.solid",false);
    bool dirtCapVisible = _config->getBool("dirt.capVisible",true);
    bool dirtCapSolid   = _config->getBool("dirt.capSolid",false);

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
    double capHalfHeight = _config->getDouble("dirt.capHalfHeight");
    double capBottomR    = _config->getDouble("dirt.capBottomRadius");
    double capTopR       = _config->getDouble("dirt.capTopRadius");

    double dsz0          = _config->getDouble("toyDS.z0");

    GeomHandle<Beamline> beamg;
    double solenoidOffset = beamg->solenoidOffset();

    G4ThreeVector dirtCapOffset( -solenoidOffset, ySurface+capHalfHeight, 
                                 dsz0+VolumeInfo::getMu2eOriginInWorld().z());

    AntiLeakRegistry& reg = edm::Service<G4Helper>()->antiLeakRegistry();
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
