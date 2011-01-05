//
// Free function to create  DS. (Detector Solenoid)
//
// $Id: constructDS.cc,v 1.1 2011/01/05 21:04:47 genser Exp $
// $Author: genser $
// $Date: 2011/01/05 21:04:47 $
//
// Original author KLG based on Mu2eWorld constructDS
//
// Notes:
// Construct the DS. Parent volume is the air inside of the hall.
// This makes 4 volumes:
//  0 - a single volume that represents the coils+cryostats in an average way.
//  1 - DS1, a small piece of DS vacuum that surrounds TS5.
//  2 - DS2, the upstream part of the DS vacuum, that has a graded field.
//  3 - DS3, the downstream part of the DS vacuum, that may have a uniform field.

// Mu2e includes.

#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/constructDS.hh"
#include "Mu2eG4/inc/nestTubs.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"


using namespace std;

namespace mu2e {

  void constructDS( const VolumeInfo& parent, 
                    SimpleConfig const * const _config
                    ){

    // Extract information from the config file.
    TubsParams detSolCoilParams( _config->getDouble("toyDS.rIn"),
                                 _config->getDouble("toyDS.rOut"),
                                 _config->getDouble("toyDS.halfLength"));

    GeomHandle<Beamline> beamg;
    double solenoidOffset = beamg->solenoidOffset();
    double rTorus         = beamg->getTS().torusRadius();
    double rCryo          = beamg->getTS().outerRadius();
    double ts5HalfLength  = beamg->getTS().getTS5().getHalfLength();

    double dsCoilZ0          = _config->getDouble("toyDS.z0");
    double ds1HalfLength     = _config->getDouble("toyDS1.halfLength");
    double ds2HalfLength     = _config->getDouble("toyDS2.halfLength");
    double ds3HalfLength     = _config->getDouble("toyDS3.halfLength");
    double dsFrontHalfLength = _config->getDouble("toyDS.frontHalfLength");

    // All Vacuum volumes fit inside the DS coil+cryostat package.
    // DS1 surrounds ts5.
    // DS2 and DS3 extend to r=0
    TubsParams dsFrontParams( rCryo,
                              detSolCoilParams.innerRadius,
                              dsFrontHalfLength);

    TubsParams ds1VacParams( rCryo,
                             detSolCoilParams.innerRadius,
                             ds1HalfLength);

    TubsParams ds2VacParams( 0.,
                             detSolCoilParams.innerRadius,
                             ds2HalfLength);

    TubsParams ds3VacParams( 0.,
                             detSolCoilParams.innerRadius,
                             ds3HalfLength);

    // Compute positions of objects in Mu2e coordinates.
    double dsFrontZ0 = rTorus + 2.*ts5HalfLength - 2.*ds1HalfLength - dsFrontHalfLength;
    double ds1Z0     = rTorus + 2.*ts5HalfLength - ds1HalfLength;
    double ds2Z0     = rTorus + 2.*ts5HalfLength + ds2HalfLength;
    double ds3Z0     = ds2Z0  + ds2HalfLength    + ds3HalfLength;
    G4ThreeVector detSolCoilPosition(-solenoidOffset, 0., dsCoilZ0);
    G4ThreeVector    dsFrontPosition(-solenoidOffset, 0., dsFrontZ0);
    G4ThreeVector        ds1Position(-solenoidOffset, 0., ds1Z0);
    G4ThreeVector        ds2Position(-solenoidOffset, 0., ds2Z0);
    G4ThreeVector        ds3Position(-solenoidOffset, 0., ds3Z0);

    MaterialFinder materialFinder(*_config);
    G4Material* detSolCoilMaterial = materialFinder.get("toyDS.materialName");
    G4Material* vacuumMaterial     = materialFinder.get("toyDS.insideMaterialName");

    // Single volume representing the DS coils + cryostat in an average way.

    bool toyDSVisible        = _config->getBool("toyDS.visible",true);
    bool toyDSSolid          = _config->getBool("toyDS.solid",true);
    bool forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;

    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

    VolumeInfo detSolCoilInfo = nestTubs( "ToyDSCoil",
                                          detSolCoilParams,
                                          detSolCoilMaterial,
                                          0,
                                          detSolCoilPosition-_hallOriginInMu2e,
                                          parent,
                                          0,
                                          toyDSVisible,
                                          G4Color::Magenta(),
                                          toyDSSolid,
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );

    // Upstream face of the DS coils+cryo.
    VolumeInfo dsFrontInfo    = nestTubs( "ToyDSFront",
                                          dsFrontParams,
                                          detSolCoilMaterial,
                                          0,
                                          dsFrontPosition-_hallOriginInMu2e,
                                          parent,
                                          0,
                                          toyDSVisible,
                                          G4Color::Blue(),
                                          toyDSSolid,
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );


    VolumeInfo ds1VacInfo     = nestTubs( "ToyDS1Vacuum",
                                          ds1VacParams,
                                          vacuumMaterial,
                                          0,
                                          ds1Position-_hallOriginInMu2e,
                                          parent,
                                          0,
                                          toyDSVisible,
                                          G4Colour::Green(),
                                          toyDSSolid,
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );

    VolumeInfo ds2VacInfo     = nestTubs( "ToyDS2Vacuum",
                                          ds2VacParams,
                                          vacuumMaterial,
                                          0,
                                          ds2Position-_hallOriginInMu2e,
                                          parent,
                                          0,
                                          toyDSVisible,
                                          G4Colour::Yellow(),
                                          toyDSSolid,
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );

    VolumeInfo ds3VacInfo     = nestTubs( "ToyDS3Vacuum",
                                          ds3VacParams,
                                          vacuumMaterial,
                                          0,
                                          ds3Position-_hallOriginInMu2e,
                                          parent,
                                          0,
                                          toyDSVisible,
                                          G4Color::Blue(),
                                          toyDSSolid,
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );

  } // end of Mu2eWorld::constructDS;

}
