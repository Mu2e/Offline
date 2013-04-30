//
// Free function to create  DS. (Detector Solenoid)
//
// $Id: constructDS.cc,v 1.11 2013/04/30 14:56:57 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/04/30 14:56:57 $
//
// Original author KLG based on Mu2eWorld constructDS
//
// Notes:
// Construct the DS. Parent volume is the air inside of the hall.
// This makes 5 volumes:
//  0 - a single volume that represents the coils+cryostats in an average way.
//  F - the front face of the coils+cryostats, made of same material as vol. 0
//  1 - DS1, a small piece of DS vacuum that surrounds TS5.
//  2 - DS2, the upstream part of the DS vacuum, that has a graded field.
//  3 - DS3, the downstream part of the DS vacuum, that may have a uniform field.
//      Length of DS3 is coupled to the length of the MBS.

// Mu2e includes.
#include "BeamlineGeom/inc/Beamline.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeomPrimitives/inc/PolyconsParams.hh"
#include "MBSGeom/inc/MBS.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/constructDS.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Polycone.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;

namespace mu2e {

  void constructDS( const VolumeInfo& parent,
                    SimpleConfig const & _config
                    ){
    
    int const verbosityLevel = _config.getInt("toyDS.verbosityLevel",0);
    
    GeomHandle<DetectorSolenoid> ds;
    TubsParams detSolCoilParams( ds->rIn(), ds->rOut(), ds->halfLength() );
    
    // All Vacuum volumes fit inside the DS coil+cryostat package.
    // - DS1 surrounds ts5.
    // - DS2 and DS3 extend to r=0    
    // - DS3 is defined as a G4Polycone object later on
    GeomHandle<Beamline> beamg;
    TubsParams dsFrontParams   ( beamg->getTS().outerRadius(), detSolCoilParams.innerRadius(), ds->frontHalfLength() );
    TubsParams ds1VacParams    ( beamg->getTS().outerRadius(), detSolCoilParams.innerRadius(), ds->halfLengthDs1()   );
    TubsParams ds2VacParams    ( 0.                          , detSolCoilParams.innerRadius(), ds->halfLengthDs2()   );

    // Compute/set positions of vacuum volumes in Mu2e coordinates.
    // - DS position is fixed by TS torus radius, and half lengths of 
    //   front face, DS1, and TS5
    double dsFrontZ0 = ds->position().z() - ds->halfLength() + ds->frontHalfLength();
    double ds1Z0     = dsFrontZ0 + ds->frontHalfLength() + ds->halfLengthDs1();
    double ds2Z0     = ds->zLocDs23Split() - ds->halfLengthDs2();

    G4ThreeVector detSolCoilPosition( ds->position().x(), 0., ds->position().z());
    G4ThreeVector    dsFrontPosition( ds->position().x(), 0., dsFrontZ0);
    G4ThreeVector        ds1Position( ds->position().x(), 0., ds1Z0);
    G4ThreeVector        ds2Position( ds->position().x(), 0., ds2Z0);

    // Check materials
    G4Material* detSolCoilMaterial = findMaterialOrThrow( ds->material() );
    G4Material* vacuumMaterial     = findMaterialOrThrow( ds->insideMaterial() );

    // Single volume representing the DS coils + cryostat in an average way.

    bool toyDSVisible        = _config.getBool("toyDS.visible",true);
    bool toyDSSolid          = _config.getBool("toyDS.solid",true);
    bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
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

    // Polycone geometry allows for MBS to extend beyond solenoid
    // physical boundaries
    
    GeomHandle<MBS> mbs;

    // Define polycone parameters
    vector<double> tmp_zPlanesDs3;
    tmp_zPlanesDs3.push_back( ds2Z0 + ds->halfLengthDs2() );
    tmp_zPlanesDs3.push_back( ds->coilZMax()              );
    tmp_zPlanesDs3.push_back( ds->coilZMax()              );
    tmp_zPlanesDs3.push_back( mbs->getEnvelopeZmax()      );

    vector<double> tmp_rOuterDs3;
    tmp_rOuterDs3.push_back( detSolCoilParams.innerRadius() );
    tmp_rOuterDs3.push_back( detSolCoilParams.innerRadius() );
    tmp_rOuterDs3.push_back( mbs->getEnvelopeRmax()         );
    tmp_rOuterDs3.push_back( mbs->getEnvelopeRmax()         );

    assert( tmp_zPlanesDs3.size() == tmp_rOuterDs3.size() );

    vector<double> tmp_rInnerDs3 ( tmp_rOuterDs3.size(), 0. );


    PolyconsParams ds3PolyParams( tmp_zPlanesDs3,
                                  tmp_rInnerDs3,
                                  tmp_rOuterDs3 );
    
    CLHEP::Hep3Vector ds3positionInMu2e( ds->position().x(), ds->position().y(), 0.);
    
    nestPolycone( "ToyDS3Vacuum",
                  ds3PolyParams,
                  vacuumMaterial,
                  0,
                  ds3positionInMu2e - parent.centerInMu2e(),
                  parent,
                  0,
                  toyDSVisible,
                  G4Colour::Yellow(),
                  toyDSSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

  } // end of Mu2eWorld::constructDS;

}
