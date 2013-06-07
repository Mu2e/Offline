//
// Free function to create DS. (Detector Solenoid)
//
// $Id: constructDS.cc,v 1.13 2013/06/07 17:59:56 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/06/07 17:59:56 $
//
// Original author KLG based on Mu2eWorld constructDS
//

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
    
    // Flags
    int const verbosityLevel = _config.getInt("ds.verbosityLevel",0);
    bool dsVisible           = _config.getBool("ds.visible",true);
    bool dsSolid             = _config.getBool("ds.solid",true);
    bool dsCoilVisible       = _config.getBool("dsCoil.visible",true);
    bool dsCoilSolid         = _config.getBool("dsCoil.solid",true);
    bool dsSupportVisible    = _config.getBool("dsSupport.visible",true);
    bool dsSupportSolid      = _config.getBool("dsSupport.solid",true);
    bool dsShieldVisible     = _config.getBool("dsShield.visible",true);
    bool dsShieldSolid       = _config.getBool("dsShield.solid",true);
    bool dsVacuumVisible     = _config.getBool("dsVacuum.visible",true);
    bool dsVacuumSolid       = _config.getBool("dsVacuum.solid",true);
    bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;

    // Fetch parent (hall) position
    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

    // Fetch DS geom. object
    GeomHandle<DetectorSolenoid> ds;
    CLHEP::Hep3Vector const & dsP ( ds->position() );

    // DS cryostat
    G4Material*   dsCryoMaterial = findMaterialOrThrow( ds->material() );

    // - inner cryo shell
    G4ThreeVector dsInnerCryoPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsInnerCryoParams  ( ds->rIn1(), ds->rIn2(), ds->halfLength() );
    nestTubs( "DSInnerCryoShell",
              dsInnerCryoParams,
              dsCryoMaterial,
              0,
              dsInnerCryoPosition-_hallOriginInMu2e,
              parent,
              0,
              dsVisible,
              G4Color::Magenta(),
              dsSolid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    // - outer cryo shell
    G4ThreeVector dsOuterCryoPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsOuterCryoParams  ( ds->rOut1(), ds->rOut2(), ds->halfLength() );
    nestTubs( "DSOuterCryoShell",
              dsOuterCryoParams,
              dsCryoMaterial,
              0,
              dsOuterCryoPosition-_hallOriginInMu2e,
              parent,
              0,
              dsVisible,
              G4Color::Magenta(),
              dsSolid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    // - end walls
    TubsParams    dsEndWallParams    ( ds->rIn2(), ds->rOut1(), ds->endWallHalfLength() );
    G4ThreeVector dsUpEndWallPosition( dsP.x(), dsP.y(), 
                                       dsP.z() - ds->halfLength() + ds->endWallHalfLength());
    nestTubs( "DSUpEndWallShell",
              dsEndWallParams,
              dsCryoMaterial,
              0,
              dsUpEndWallPosition-_hallOriginInMu2e,
              parent,
              0,
              dsVisible,
              G4Color::Magenta(),
              dsSolid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    G4ThreeVector dsDownEndWallPosition( dsP.x(), dsP.y(), 
                                         dsP.z() + ds->halfLength() - ds->endWallHalfLength());
    nestTubs( "DSDownEndWallShell",
              dsEndWallParams,
              dsCryoMaterial,
              0,
              dsDownEndWallPosition-_hallOriginInMu2e,
              parent,
              0,
              dsVisible,
              G4Color::Magenta(),
              dsSolid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    // - upstream face
    GeomHandle<Beamline> beamg;
    double dsFrontZ0 = dsP.z() - ds->halfLength() + ds->frontHalfLength();
    TubsParams    dsFrontParams  ( beamg->getTS().outerRadius(), dsInnerCryoParams.innerRadius(), ds->frontHalfLength() );
    G4ThreeVector dsFrontPosition( dsP.x(), dsP.y(), dsFrontZ0);

    nestTubs( "DSFront",
              dsFrontParams,
              dsCryoMaterial,
              0,
              dsFrontPosition-_hallOriginInMu2e,
              parent,
              0,
              dsVisible,
              G4Color::Blue(),
              dsSolid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    // DS thermal shield
    G4Material*   dsShieldMaterial = findMaterialOrThrow( ds->shield_material() );

    // - inner shield shell
    G4ThreeVector dsInnerShieldPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsInnerShieldParams  ( ds->shield_rIn1(), ds->shield_rIn2(), ds->shield_halfLength() );
    nestTubs( "DSInnerShieldShell",
              dsInnerShieldParams,
              dsShieldMaterial,
              0,
              dsInnerShieldPosition-_hallOriginInMu2e,
              parent,
              0,
              dsShieldVisible,
              G4Color::Cyan(),
              dsShieldSolid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    // - outer shield shell
    G4ThreeVector dsOuterShieldPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsOuterShieldParams  ( ds->shield_rOut1(), ds->shield_rOut2(), ds->shield_halfLength() );
    nestTubs( "DSOuterShieldShell",
              dsOuterShieldParams,
              dsShieldMaterial,
              0,
              dsOuterShieldPosition-_hallOriginInMu2e,
              parent,
              0,
              dsShieldVisible,
              G4Color::Cyan(),
              dsShieldSolid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    // SKIPPING END WALLS OF THERMAL SHIELD
    // DUE TO INCONSISTENCY IN DRAWINGS!

    // DS coils
    G4Material* dsCoilMaterial = findMaterialOrThrow( ds->coil_material() );
    
    for ( int i(0); i < ds->nCoils() ; i++ ) {
      TubsParams coilParams( ds->coil_rIn(), 
                             ds->coil_rOut().at(i),
                             ds->coil_zLength().at(i)*0.5 );
      G4ThreeVector coilPosition( dsP.x(), dsP.y(), ds->coil_zPosition().at(i) );

      ostringstream coilname;
      coilname << "DSCoil_" << i+1;

      nestTubs( coilname.str(),
                coilParams,
                dsCoilMaterial,
                0,
                coilPosition-_hallOriginInMu2e,
                parent,
                0,
                dsCoilVisible,
                G4Color::Green(),
                dsCoilSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
    }

    // DS coils support system
    G4Material*   dsSupportMaterial = findMaterialOrThrow( ds->support_material() );
    G4ThreeVector dsSupportPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsSupportParams  ( ds->support_rIn(), ds->support_rOut(), ds->support_halfLength() );
    nestTubs( "DSCoilSupport",
              dsSupportParams,
              dsSupportMaterial,
              0,
              dsSupportPosition-_hallOriginInMu2e,
              parent,
              0,
              dsSupportVisible,
              G4Color::Blue(),
              dsSupportSolid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    // DS vacuum volumes
    G4Material* vacuumMaterial = findMaterialOrThrow( ds->vacuumMaterial() );
    TubsParams ds1VacParams    ( beamg->getTS().outerRadius(), ds->rIn1(), ds->vac_halfLengthDs1()   );
    TubsParams ds2VacParams    ( 0.                          , ds->rIn1(), ds->vac_halfLengthDs2()   );

    // Compute/set positions of vacuum volumes in Mu2e coordinates.
    // - DS position is fixed by TS torus radius, and half lengths of 
    //   front face, DS1, and TS5
    double ds1Z0     = dsFrontZ0 + ds->frontHalfLength() + ds->vac_halfLengthDs1();
    double ds2Z0     = ds->vac_zLocDs23Split() - ds->vac_halfLengthDs2();

    G4ThreeVector ds1Position( dsP.x(), dsP.y(), ds1Z0 );
    G4ThreeVector ds2Position( dsP.x(), dsP.y(), ds2Z0 );

    nestTubs( "DS1Vacuum",
              ds1VacParams,
              vacuumMaterial,
              0,
              ds1Position-_hallOriginInMu2e,
              parent,
              0,
              dsVacuumVisible,
              G4Colour::Green(),
              dsVacuumSolid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    nestTubs( "DS2Vacuum",
              ds2VacParams,
              vacuumMaterial,
              0,
              ds2Position-_hallOriginInMu2e,
              parent,
              0,
              dsVacuumVisible,
              G4Colour::Yellow(),
              dsVacuumSolid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    // Polycone geometry allows for MBS to extend beyond solenoid
    // physical boundaries
    
    GeomHandle<MBS> mbs;

    // Define polycone parameters
    vector<double> tmp_zPlanesDs3 { ds->vac_zLocDs23Split(), ds->cryoZMax(), ds->cryoZMax(),         mbs->getEnvelopeZmax() };
    vector<double> tmp_rOuterDs3  { ds->rIn1(),              ds->rIn1(),     mbs->getEnvelopeRmax(), mbs->getEnvelopeRmax() };
    assert( tmp_zPlanesDs3.size() == tmp_rOuterDs3.size() );

    vector<double> tmp_rInnerDs3 ( tmp_rOuterDs3.size(), 0. );

    PolyconsParams    ds3PolyParams    ( tmp_zPlanesDs3, tmp_rInnerDs3, tmp_rOuterDs3 );
    CLHEP::Hep3Vector ds3positionInMu2e( dsP.x(), dsP.y(), 0.);
    
    nestPolycone( "DS3Vacuum",
                  ds3PolyParams,
                  vacuumMaterial,
                  0,
                  ds3positionInMu2e - parent.centerInMu2e(),
                  parent,
                  0,
                  dsVacuumVisible,
                  G4Colour::Yellow(),
                  dsVacuumSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

  } // end of Mu2eWorld::constructDS;

}
