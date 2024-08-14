//
// Free function to create DS. (Detector Solenoid)
//
//
// Original author KLG based on Mu2eWorld constructDS
//

// Mu2e includes.
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/BeamlineGeom/inc/StraightSection.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoidShielding.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeomPrimitives/inc/PolyconsParams.hh"
#include "Offline/MBSGeom/inc/MBS.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/constructDS.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestPolycone.hh"
#include "Offline/Mu2eG4/inc/nestExtrudedSolid.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/Mu2eG4/inc/checkForOverlaps.hh"

// G4 includes
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4ExtrudedSolid.hh"
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4SDManager.hh"
// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
// art includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

using namespace std;

namespace mu2e {

  void constructDS( const VolumeInfo& parent,
                    const SimpleConfig& _config
                    ){
    MaterialFinder materialFinder(_config);

    // Load flags
    int const verbosityLevel = _config.getInt("ds.verbosityLevel",0);
    bool const inGaragePosition = _config.getBool("inGaragePosition",false);
    double zOffGarage = (inGaragePosition) ? _config.getDouble("garage.zOffset") : 0.;
    bool const OPA_IPA_ST_Extracted = (inGaragePosition) ? _config.getBool("garage.extractOPA_IPA_ST") : false;
    //insert Z offset for the extracted position
    CLHEP::Hep3Vector relPosExtracted(0.,0., zOffGarage);

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "ds"         , "ds"           );
    geomOptions->loadEntry( _config, "dsCoil"     , "ds.coil"      );
    geomOptions->loadEntry( _config, "dsRing"     , "ds.ring"      );
    geomOptions->loadEntry( _config, "dsSpacer"   , "ds.spacer"    );
    geomOptions->loadEntry( _config, "dsSupport"  , "ds.support"   );
    geomOptions->loadEntry( _config, "dsThShield" , "ds.thShield"  );
    geomOptions->loadEntry( _config, "dsVacuum"   , "ds.vacuum"    );
    geomOptions->loadEntry( _config, "dsShielding", "ds.shielding" );

    const bool isDSVisible         = geomOptions->isVisible("ds");
    const bool isDSSolid           = geomOptions->isSolid("ds");
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("ds");
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck("ds");
    const bool placePV             = geomOptions->placePV("ds");

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

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
              G4Color::Magenta(),
              "ds"
              );
    // ***
    // Lining for inner cryo shell - for test of shielding
    // added 11 June 2017
    // ***
    if ( ds->hasInnerLining() ) { // only if specifically enabled
      nestTubs( "DSInnerCryoLining",
                TubsParams( ds->rIn1()-ds->innerLiningThickness(),ds->rIn1(),ds->halfLength()-2.*ds->endWallHalfLength()),
                findMaterialOrThrow(ds->innerLiningMaterial()),
                0,
                dsInnerCryoPosition-_hallOriginInMu2e,
                parent,
                0,
                G4Color::Magenta(),
                "ds"
                );
    }

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
              G4Color::Magenta(),
              "ds"
              );

    // DNB (Lou) May 2017.  Making a vacuum volume inside the cryostat
    G4Tubs * dsCryoVacTub = new G4Tubs( "DSCryoVacuumTube",
                                        ds->rIn2(), ds->rOut1(),
                                        ds->halfLength() - 2.*ds->endWallHalfLength(),
                                        0.0, 360.0*CLHEP::degree);

    CLHEP::Hep3Vector dsCryoVacLocationInMu2e( dsP.x(), dsP.y(), dsP.z() );
    VolumeInfo dsCryoVacMother( "DSCryoVacuumRegion",
                                dsCryoVacLocationInMu2e - _hallOriginInMu2e,
                                parent.centerInWorld );
    dsCryoVacMother.solid = dsCryoVacTub;

    finishNesting ( dsCryoVacMother,
                    findMaterialOrThrow("DSVacuum"),
                    0, dsCryoVacMother.centerInParent,
                    parent.logical, 0, G4Colour::White(),
                    "ds" );

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
              G4Color::Magenta(),
              "ds"
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
              G4Color::Magenta(),
              "ds"
              );

    // - upstream face
    GeomHandle<Beamline> beamg;
    auto ts(&beamg->getTS());
    const StraightSection * ts5out = ts->getTSCryo<StraightSection>( TransportSolenoid::TSRegion::TS5,
                                                                                TransportSolenoid::TSRadialPart::OUT );

    double interconnectHL = _config.getDouble("ds.interconnect.halfLength", -1.);
    double dsFrontZ0 = dsP.z() - ds->halfLength() + ds->frontHalfLength();
    double smallOffset = 1.e-3; //create small buffers to prevent overlaps on both ends of the interconnect
    if(interconnectHL < 0.) {
      TubsParams    dsFrontParams  ( ts5out->rOut(), dsInnerCryoParams.innerRadius(), ds->frontHalfLength() );
      G4ThreeVector dsFrontPosition( dsP.x(), dsP.y(), dsFrontZ0);

      nestTubs( "DSFront",
                dsFrontParams,
                dsCryoMaterial,
                0,
                dsFrontPosition-_hallOriginInMu2e,
                parent,
                0,
                G4Color::Blue(),
                "ds"
                );
    } else {
      double interconnectRIn  = _config.getDouble("ds.interconnect.rIn");
      double interconnectROut = _config.getDouble("ds.interconnect.rOut");
      nestTubs( "TSDSInterConnect",
                TubsParams(interconnectRIn, interconnectROut, interconnectHL-smallOffset),
                dsCryoMaterial,
                0,
                dsP-_hallOriginInMu2e-CLHEP::Hep3Vector(0.,0.,ds->halfLength()+interconnectHL+smallOffset),
                parent,
                0,
                G4Color::Blue(),
                "ds"
                );
    }

    // DS thermal shield
    G4Material*   dsThShieldMaterial = findMaterialOrThrow( ds->shield_material() );

    // - inner shield shell
    G4ThreeVector dsInnerThShieldPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsInnerThShieldParams  ( ds->shield_rIn1(), ds->shield_rIn2(), ds->shield_halfLength() );
    nestTubs( "DSInnerThShieldShell",
              dsInnerThShieldParams,
              dsThShieldMaterial,
              0,
              dsInnerThShieldPosition-dsCryoVacLocationInMu2e,
              dsCryoVacMother,
              0,
              G4Color::Cyan(),
              "dsThShield"
              );

    // - outer shield shell
    G4ThreeVector dsOuterThShieldPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsOuterThShieldParams  ( ds->shield_rOut1(), ds->shield_rOut2(), ds->shield_halfLength() );
    nestTubs( "DSOuterThShieldShell",
              dsOuterThShieldParams,
              dsThShieldMaterial,
              0,
              dsOuterThShieldPosition-dsCryoVacLocationInMu2e,
              dsCryoVacMother,
              0,
              G4Color::Cyan(),
              "dsThShield"
              );

    // SKIPPING END WALLS OF THERMAL SHIELD
    // DUE TO INCONSISTENCY IN DRAWINGS!

    // DS coils
    G4Material* dsCoilMaterial;
    for ( int i(0); i < ds->nCoils() ; i++ ) {
      TubsParams coilParams( ds->coil_rIn(),
                             ds->coil_rOut().at(i),
                             ds->coil_zLength().at(i)*0.5 );
      if ( ds->coilVersion() == 1 ) {
        dsCoilMaterial = findMaterialOrThrow( ds->coil_material() );
      } else {
        dsCoilMaterial = findMaterialOrThrow( ds->coil_materials().at(i) );
      }

      G4ThreeVector coilPosition( dsP.x(), dsP.y(),
                                  ds->coil_zPosition().at(i) + coilParams.zHalfLength());

      std::ostringstream coilname;
      coilname << "DSCoil_" << i+1;

      nestTubs( coilname.str(),
                coilParams,
                dsCoilMaterial,
                0,
                coilPosition-dsCryoVacLocationInMu2e,
                dsCryoVacMother,
                0,
                G4Color::Green(),
                "dsCoil"
                );
    }

    // DS coil spacers
    if ( ds->coilVersion() > 1 ) {  // spacers added to version 2 and up
      G4Material* dsSpacerMaterial = findMaterialOrThrow( ds->spacer_material() );
      for ( int i(0); i < ds->nSpacers() ; i++ ) {
        TubsParams spacerParams( ds->spacer_rIn(),
                                 ds->spacer_rOut().at(i),
                                 ds->spacer_zLength().at(i)*0.5 );

        G4ThreeVector spacerPosition( dsP.x(), dsP.y(),
                                      ds->spacer_zPosition().at(i)
                                      + spacerParams.zHalfLength());

        ostringstream spacername;
        spacername << "DSSpacer_" << i+1;

        nestTubs( spacername.str(),
                  spacerParams,
                  dsSpacerMaterial,
                  0,
                  spacerPosition-dsCryoVacLocationInMu2e,
                  dsCryoVacMother,
                  0,
                  G4Color::Green(),
                  "dsSpacer"
                  );
      }
    }


    // DS coils support system
    G4Material*   dsSupportMaterial = findMaterialOrThrow( ds->support_material() );
    G4ThreeVector dsSupportPosition( dsP.x(), dsP.y(), dsP.z());
    TubsParams    dsSupportParams  ( ds->support_rIn(), ds->support_rOut(), ds->support_halfLength() );
    nestTubs( "DSCoilSupport",
              dsSupportParams,
              dsSupportMaterial,
              0,
              dsSupportPosition-dsCryoVacLocationInMu2e,
              dsCryoVacMother,
              0,
              G4Color::Blue(),
              "dsSupport"
              );



    // Build DS Rings.  Added by David Norvil Brown, May 2015
    double rirs = ds->rInRingSide();
    double rors = ds->rOutRingSide();
    double trs = ds->thickRingSide();
    double rir = ds->rInRing();
    double ror = ds->rOutRing();
    double lr = ds->lengthRing();
    G4Material* ringMaterial = findMaterialOrThrow(ds->RingMaterial());
    std::vector<double> xr = ds->xRing();
    std::vector<double> yr = ds->yRing();
    std::vector<double> zr = ds->zRing();

    for ( unsigned int iRing = 0; iRing < xr.size(); iRing++ ) {

      // Let's build a mother volume first
      std::ostringstream ringMotherName;
      ringMotherName << "DSRingMother" << iRing;

      double motherx = xr[iRing];
      double mothery = yr[iRing];
      double motherz = zr[iRing];

      VolumeInfo motherVol = nestTubs( ringMotherName.str(),
                                       TubsParams( rirs, rors, trs + lr/2.0 ),
                                       findMaterialOrThrow("G4_AIR"),
                                       nullptr,
                                       CLHEP::Hep3Vector(motherx,mothery,motherz) - _hallOriginInMu2e,
                                       parent, 0, G4Color::Blue(),
                                       "dsRing" );

      std::ostringstream leftName;
      leftName << "DSleftSideRing" << iRing;


      nestTubs( leftName.str(),
                TubsParams( rirs, rors, trs/2.0 ),
                ringMaterial,
                nullptr,
                CLHEP::Hep3Vector(0.0,0.0,-lr/2.0-trs/2.0),
                motherVol,
                0,
                G4Color::Blue(),
                "dsRing"
                );

      std::ostringstream centerName;
      centerName << "DScenterRing" << iRing;

      nestTubs( centerName.str(),
                TubsParams( rir, ror, lr/2.0 ),
                ringMaterial,
                nullptr,
                CLHEP::Hep3Vector(0.0,0.0,0.0),
                motherVol,
                0,
                G4Color::Blue(),
                "dsRing"
                );

      std::ostringstream rightName;
      rightName << "DSrightSideRing" << iRing;

      nestTubs( rightName.str(),
                TubsParams( rirs, rors, trs/2.0 ),
                ringMaterial,
                nullptr,
                CLHEP::Hep3Vector(0.0,0.0,lr/2.0+trs/2.0),
                motherVol,
                0,
                G4Color::Blue(),
                "dsRing"
                );

    } // finished inserting Rings


    // DS vacuum volumes
    G4Material* vacuumMaterial = findMaterialOrThrow( ds->vacuumMaterial() );
    double extendDS1Vacuum = 0.;
    if(interconnectHL >= 0.) extendDS1Vacuum = interconnectHL + ds->frontHalfLength() - smallOffset;
    TubsParams ds1VacParams    ( ts5out->rOut(), ds->rIn1(), ds->vac_halfLengthDs1() + extendDS1Vacuum  );
    TubsParams ds2VacParams    ( 0.            , ds->rIn1(), ds->vac_halfLengthDs2()   );

    // Compute/set positions of vacuum volumes in Mu2e coordinates.
    // - DS position is fixed by TS torus radius, and half lengths of
    //   front face, DS1, and TS5
    double ds1Z0     = dsFrontZ0 + ds->frontHalfLength() + ds->vac_halfLengthDs1() - extendDS1Vacuum;
    double ds2Z0     = ds->vac_zLocDs23Split() - ds->vac_halfLengthDs2();
    //    double ds2HalfLength     = _config.getDouble("ds2.halfLength");

    if ( verbosityLevel > 0 ) {
      G4cout << __func__ << " DS2 vacuum extent: "
           << " [ " << ds2Z0 - ds->vac_halfLengthDs2() << " , "
           << ds->vac_zLocDs23Split() << " ] " << G4endl;
    }

    G4ThreeVector ds1Position( dsP.x(), dsP.y(), ds1Z0 );
    G4ThreeVector ds2Position( dsP.x(), dsP.y(), ds2Z0 );

    nestTubs( "DS1Vacuum",
              ds1VacParams,
              vacuumMaterial,
              0,
              ds1Position-_hallOriginInMu2e,
              parent,
              0,
              G4Colour::Green(),
              "dsVacuum"
              );

    VolumeInfo ds2VacInfo =
      nestTubs( "DS2Vacuum",
                ds2VacParams,
                vacuumMaterial,
                0,
                ds2Position-_hallOriginInMu2e,
                parent,
                0,
                G4Colour::Yellow(),
                "dsVacuum"
                );

    VolumeInfo tmpDS;
    //create volume for detector elements in the extracted position
    if(inGaragePosition && OPA_IPA_ST_Extracted) {
      G4Material*  airMaterial = findMaterialOrThrow( _config.getString("hall.insideMaterialName","G4_AIR") );

      VolumeInfo dsShieldParent = nestTubs( "garageFakeDS2Vacuum",
                                            ds2VacParams,
                                            airMaterial,
                                            0,
                                            ds2Position - _hallOriginInMu2e + relPosExtracted,
                                            parent,
                                            0,
                                            G4Colour::Yellow(),
                                            "dsVacuum"
                                            );
    }
    // Polycone geometry allows for MBS to extend beyond solenoid
    // physical boundaries

    GeomHandle<DetectorSolenoidShielding> dss;

    // Define polycone parameters
    vector<double> tmp_zPlanesDs3 = {
      ds->vac_zLocDs23Split(),
      dss->getVPSPmain()->zBegin(),
      dss->getVPSPCryoSeal()->zBegin(),
      dss->getVPSPCryoSeal()->zEnd(),
      dss->getVPSPCryoSeal()->zEnd(),
      dss->getVPSPendSeal()->zBegin(),
      dss->getVPSPendSeal()->zBegin(),
      dss->getVPSPendFlange()->zEnd(),
      dss->getVPSPendFlange()->zEnd(),
      dss->getIFBmain()->zEnd(),
      dss->getIFBendPlug()->zBegin(),
      dss->getIFBendPlug()->zEnd()
    };

    vector<double> tmp_rOuterDs3  = {
      ds->rIn1(),
      dss->getVPSPmain()->innerRadius(),
      dss->getVPSPCryoSeal()->outerRadius(),
      dss->getVPSPCryoSeal()->outerRadius(),
      dss->getVPSPmain()->outerRadius(),
      dss->getVPSPmain()->outerRadius(),
      dss->getVPSPendSeal()->outerRadius(),
      dss->getVPSPendFlange()->outerRadius(),
      dss->getVPSPendFlange()->innerRadius(),
      dss->getIFBmain()->outerRadius(),
      dss->getIFBendPlug()->outerRadius(),
      dss->getIFBendPlug()->outerRadius()
    };

    assert( tmp_zPlanesDs3.size() == tmp_rOuterDs3.size() );

    vector<double> tmp_rInnerDs3 ( tmp_rOuterDs3.size(), 0. );

    PolyconsParams    ds3PolyParams    ( tmp_zPlanesDs3, tmp_rInnerDs3, tmp_rOuterDs3 );
    CLHEP::Hep3Vector ds3positionInMu2e( dsP.x(), dsP.y(), 0.);

    VolumeInfo dsShieldParent = nestPolycone( "DS3Vacuum",
                                              ds3PolyParams,
                                              vacuumMaterial,
                                              0,
                                              ds3positionInMu2e - parent.centerInMu2e(),
                                              parent,
                                              0,
                                              G4Colour::Yellow(),
                                              "dsVacuum"
                                              );
    if ( inGaragePosition ) {
      G4Material*  airMaterial = findMaterialOrThrow( _config.getString("hall.insideMaterialName","G4_AIR") );

      tmpDS = nestPolycone( "garageFakeDS3Vacuum",
                            ds3PolyParams,
                            airMaterial,
                            0,
                            ds3positionInMu2e - parent.centerInMu2e() + relPosExtracted,
                            parent,
                            0,
                            G4Colour::Yellow(),
                            "dsVacuum"
                            );
    } else {
      tmpDS = dsShieldParent;
    }
    VolumeInfo & dsShieldPointer = tmpDS;


    // Construct shielding downstream of DS
    for ( const auto & shield : dss->getTubes() ) {

      nestTubs( shield->getName(),
                shield->getTubsParams(),
                findMaterialOrThrow(shield->materialName()),
                0,
                shield->originInMu2e()-dsShieldParent.centerInMu2e(),
                dsShieldParent,
                0,
                G4Colour::Blue(),
                "dsShielding"
                );

      if (verbosityLevel > 0) {
        cout << __func__ << ": " << shield->getName() << " Params = " << shield->getTubsParams() << ", location = " << shield->originInMu2e() << endl;
      }
    }

    // ************ End Shielding, begin Rails ************
    // D. No. Brown
    // Add rails in DS2Vacuum

    //    std::vector<CLHEP::Hep2Vector> railOutline;
    std::vector<double> uRailOutline = ds->uOutlineRail();
    std::vector<double> vRailOutline = ds->vOutlineRail();

    CLHEP::HepRotation* nRailRotat = nullptr;
    CLHEP::HepRotation* sRailRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
    sRailRotat->rotateY(180.0*CLHEP::degree);

    VolumeInfo RailN2 = nestExtrudedSolid
                     ( "NorthRailDS2", ds->lengthRail2()/2.0*CLHEP::mm,
                       uRailOutline, vRailOutline,
                       findMaterialOrThrow(ds->RailMaterial()),
                       nRailRotat, ds->n2RailCenter(),
                       ds2VacInfo.logical, 0, isDSVisible,
                       G4Colour::Blue(), isDSSolid,
                       forceAuxEdgeVisible, placePV, doSurfaceCheck );

    VolumeInfo RailS2 = nestExtrudedSolid
                     ( "SouthRailDS2", ds->lengthRail2()/2.0*CLHEP::mm,
                       uRailOutline, vRailOutline,
                       findMaterialOrThrow(ds->RailMaterial()),
                       sRailRotat, ds->s2RailCenter(),
                       ds2VacInfo.logical, 0, isDSVisible,
                       G4Colour::Blue(), isDSSolid,
                       forceAuxEdgeVisible, placePV, doSurfaceCheck );

    // And now in DS3Vacuum

     VolumeInfo RailN3 = nestExtrudedSolid
                      ( "NorthRailDS3", ds->lengthRail3()/2.0*CLHEP::mm,
                       uRailOutline, vRailOutline,
                       findMaterialOrThrow(ds->RailMaterial()),
                       nRailRotat, ds->n3RailCenter(),
                       dsShieldParent, 0, isDSVisible,
                       G4Colour::Blue(), isDSSolid,
                       forceAuxEdgeVisible, placePV, doSurfaceCheck );

     VolumeInfo RailS3 = nestExtrudedSolid
                      ( "SouthRailDS3", ds->lengthRail3()/2.0*CLHEP::mm,
                       uRailOutline, vRailOutline,
                       findMaterialOrThrow(ds->RailMaterial()),
                       sRailRotat, ds->s3RailCenter(),
                       dsShieldParent, 0, isDSVisible,
                       G4Colour::Blue(), isDSSolid,
                       forceAuxEdgeVisible, placePV, doSurfaceCheck );

     // Now put bearing blocks on rails
     // D. No. Brown, Jan 2016
     // First in DS2Vacuum region
     std::vector<double> uBBlockOutline = ds->uOutlineBBlock();
     std::vector<double> vBBlockOutline = ds->vOutlineBBlock();
     CLHEP::HepRotation* BBRotat = nullptr;
     double lenBB2 = ds->lengthBBlock2()/2.0*CLHEP::mm;
     double lenBB3 = ds->lengthBBlock3()/2.0*CLHEP::mm;
     std::vector<CLHEP::Hep3Vector> BBCenters2 = ds->BBlockCenters2();
     std::vector<CLHEP::Hep3Vector> BBCenters3 = ds->BBlockCenters3();
     int cScheme = ds->couplerScheme();

     // First in DS2Vacuum region
     int nB2 = BBCenters2.size();
     CLHEP::Hep3Vector DS2Offset(0,0,ds2Z0);
     int coupleCounter = 0;
     double widCoupler = ds->widthCoupler();
     double hCoupler = ds->heightCoupler();
     double yCoupler = ds->yCenterCoupler();
     for ( int iB2 = 0; iB2 < nB2; iB2++ ) {
       std::stringstream sstm;
       sstm << "BearingBlock_DS2_" << iB2+1;
       VolumeInfo BBlock2 = nestExtrudedSolid
         ( sstm.str().c_str(), lenBB2,
           uBBlockOutline, vBBlockOutline,
           findMaterialOrThrow(ds->BBlockMaterial()),
           BBRotat, BBCenters2[iB2] - DS2Offset,
           ds2VacInfo.logical,
           0, isDSVisible,
           G4Colour::Blue(), isDSSolid,
           forceAuxEdgeVisible, placePV, doSurfaceCheck );
       // Now add Couplers
       if (iB2 < nB2 - 2 && (cScheme == 0 || (cScheme == 1 && BBCenters2[iB2].x() > 0.0) || (cScheme == 2 && BBCenters2[iB2].x() < 0) ) ) {
         coupleCounter++;
         std::stringstream couplerName;
         couplerName << "Coupler_DS2_" << coupleCounter;
         double lenCoupler = BBCenters2[iB2+2].z() - BBCenters2[iB2].z()
           - 2.*lenBB2 - 0.2; // The 0.2 is to avoid accidental overlaps.
         // The 2.* is because lenBB2 is actually halfLength
         CLHEP::Hep3Vector cenCoupler( (BBCenters2[iB2] + BBCenters2[iB2+2] )*0.5 + CLHEP::Hep3Vector(0.0,yCoupler,0.0));
         std::vector<double> halfDims = { widCoupler/2.0,
                                          hCoupler/2.0,
                                          lenCoupler/2.0 };
         nestBox( couplerName.str(),
                  halfDims,
                  findMaterialOrThrow(ds->BBlockMaterial()),
                  BBRotat,
                  cenCoupler - DS2Offset,
                  ds2VacInfo.logical,
                  0, isDSVisible,
                  G4Colour::Blue(), isDSSolid,
                  forceAuxEdgeVisible, placePV, doSurfaceCheck );

       } // end of if for adding coupler if not last bearing block
     }

     // Now in DS3Vacuum region
     int nB3 = BBCenters3.size();
     coupleCounter = 0;
     for ( int iB3 = 0; iB3 < nB3; iB3++ ) {
       std::stringstream sstm;
       sstm << "BearingBlock_DS3_" << iB3+1;

       VolumeInfo BBlock3 = nestExtrudedSolid
         ( sstm.str().c_str(), lenBB3,
           uBBlockOutline, vBBlockOutline,
           findMaterialOrThrow(ds->BBlockMaterial()),
           BBRotat, BBCenters3[iB3],
           dsShieldParent, 0, isDSVisible,
           G4Colour::Blue(), isDSSolid,
           forceAuxEdgeVisible, placePV, doSurfaceCheck );
       // Now add Couplers
       if (iB3 < nB3 - 2 && (cScheme == 0 || (cScheme == 1 && BBCenters3[iB3].x() > 0.0) || (cScheme == 2 && BBCenters3[iB3].x() < 0) ) ) {
         coupleCounter++;
         std::stringstream couplerName;
         couplerName << "Coupler_DS3_" << coupleCounter;
         double lenCoupler = BBCenters3[iB3+2].z() - BBCenters3[iB3].z()
           - 2.*lenBB3 - 0.2; // The 0.2 is to avoid accidental overlaps.
         // The 2.* is because lenBB3 is actually halfLength
         CLHEP::Hep3Vector cenCoupler( (BBCenters3[iB3] + BBCenters3[iB3+2] )*0.5 + CLHEP::Hep3Vector(0.0,yCoupler,0.0));
         std::vector<double> halfDims = { widCoupler/2.0,
                                          hCoupler/2.0,
                                          lenCoupler/2.0 };
         nestBox( couplerName.str(),
                  halfDims,
                  findMaterialOrThrow(ds->BBlockMaterial()),
                  BBRotat,
                  cenCoupler,
                  dsShieldParent,
                  0, isDSVisible,
                  G4Colour::Blue(), isDSSolid,
                  forceAuxEdgeVisible, placePV, doSurfaceCheck );

       } // end of if for adding coupler if not last bearing block


     }


     //************* End of Rails and Bearing Blocks ************
     //************* Begin MBS Support block *****************
     if ( ds->hasMBSS() ) {
       std::vector<double> uOutlineMBSS = ds->uOutlineMBSS();
       std::vector<double> vOutlineMBSS = ds->vOutlineMBSS();
       VolumeInfo MBSS = nestExtrudedSolid
         ( "MBSSphericalSupport", ds->lengthMBSS()/2.0*CLHEP::mm,
           uOutlineMBSS, vOutlineMBSS,
           findMaterialOrThrow(ds->MBSSmaterial()),
           0, ds->MBSSlocation(),
           dsShieldParent, 0, isDSVisible,
           G4Colour::Blue(), isDSSolid,
           forceAuxEdgeVisible, placePV, doSurfaceCheck );
     }  // end of if ( ds->hasMBSS() )

     // End of MBS spherical shielding, begin cable runs for Cal and Tracker
     // Each is modeled as a thin wedge of a ring

     //bool cableRunSensitive = _config.getBool("ds.CableRun.sensitive",false);

     if ( ds->hasCableRunCal() ) {

       // fixme check if one should use  ds->cableRunVersion() > 1
       TubsParams  calCableRunParams  ( ds->rInCableRunCal(),
                                        ds->rOutCableRunCal(),
                                        ds->lengthCableRunCal(),
                                        ds->phi0CableRunCal()*CLHEP::degree,
                                        ds->dPhiCableRunCal()*CLHEP::degree);

       CLHEP::Hep3Vector calCableRunLoc( 0.0, 0.0, ds->zCCableRunCal() );
       VolumeInfo ccrTemp = nestTubs( "CalCableRun",
                                      calCableRunParams,
                                      findMaterialOrThrow(ds->calCableRunMaterial()),
                                      0,
                                      calCableRunLoc,
                                      dsShieldPointer,
                                      0,
                                      G4Color::Magenta(),
                                      "ds"
                                      );

       if ( ds->cableRunVersion() > 2 ) {

         // "Fibre Core"
         placeTubeCore ( "CalCableRunCore",
                         ds->rCableRunCalCoreFract(),
                         ds->rdCableRunCalCoreFract(),
                         ds->dPhiCableRunCalCoreFract(),
                         ds->materialCableRunCalCore(),
                         G4Color::Yellow(),
                         ccrTemp,
                         calCableRunParams,
                         "ds",
                         _config
                         );

       }



       if ( ds->cableRunVersion() > 1) {

         // And last but not least the connector between the top of the Cal
         // and the top of the MBS
         // Implement this as a Polycone
         std::vector<double> zs = { ds->upZC2CableRunCal() + ds->upHL2CableRunCal() + 4.0, ds->zCCableRunCal() - ds->lengthCableRunCal() };
         std::vector<double> rins = { ds->upRInCableRunCal(), ds->rInCableRunCal()};
         std::vector<double> routs= { ds->upROutCableRunCal(), ds->rOutCableRunCal()};
         PolyconsParams myPars( zs, rins, routs,
                                ds->phi0CableRunCal()*CLHEP::degree,
                                ds->dPhiCableRunCal()*CLHEP::degree );

         VolumeInfo ccrTmpF = nestPolycone ( "calCableRunFall",
                                             myPars,
                                             findMaterialOrThrow(ds->calCableRunMaterial()),
                                             0,
                                             G4ThreeVector(0,0,0),
                                             dsShieldPointer,
                                             0,
                                             G4Colour::Magenta(),
                                             "ds" );


         if ( ds->cableRunVersion() > 2 ) {

           // "Fibre Core"
           // similar to placeTubeCore below but
           // using PolyconsParams

           std::vector<double> coreInnerRadii;
           coreInnerRadii.reserve(rins.size());
           std::vector<double> coreOuterRadii;
           coreOuterRadii.reserve(rins.size());

           // loop over rins, routs

           for (  std::vector<double>::size_type i=0; i !=rins.size() ; ++i ) {

             double coreCenterRadius = routs[i] -
               ( routs[i] - rins[i] ) * ds->rCableRunCalCoreFract();
             // 0 -> rOut; 1 -> rIn; 0.5 -> (rIn + rOut)/2;

             double coreRadialHalfExtent =
               ( routs[i] - rins[i] ) * 0.5 *
               ds->rdCableRunCalCoreFract();

             coreInnerRadii.push_back( coreCenterRadius - coreRadialHalfExtent );
             coreOuterRadii.push_back( coreCenterRadius + coreRadialHalfExtent );

             if ( coreInnerRadii[i] < rins[i] ||
                  coreOuterRadii[i] > routs[i] ) {

               throw cet::exception("GEOM") << __func__
                                            << " inconsitent cable core parameters: "
                                            << i << ": "
                                            << coreInnerRadii[i]
                                            << ", "
                                            << coreOuterRadii[i]
                                            << "\n";

             }

           }

           // angular myPars are multiplied by CLHEP::degree already

           double coreCenterPhi = myPars.phi0() + 0.5 * myPars.phiTotal();
           double coreAngularHalfExtent = myPars.phiTotal() * ds->dPhiCableRunCalCoreFract() * 0.5;
           double corePhi0 = coreCenterPhi - coreAngularHalfExtent;
           double coreDeltaPhi = 2.0 * coreAngularHalfExtent;

           if ( corePhi0 < myPars.phi0() ) {

             throw cet::exception("GEOM") << __func__
                                          << " inconsitent cable core parameters: "
                                          << corePhi0
                                          << ", "
                                          << coreDeltaPhi
                                          << "\n";

           }

           PolyconsParams cableRunCoreParams (zs, coreInnerRadii, coreOuterRadii,
                                              corePhi0,
                                              coreDeltaPhi );

           VolumeInfo ccrTmpFCore = nestPolycone ( "calCableRunFallCore",
                                                   cableRunCoreParams,
                                                   findMaterialOrThrow(ds->materialCableRunCalCore()),
                                                   0,
                                                   G4ThreeVector(0,0,0),
                                                   ccrTmpF,
                                                   0,
                                                   G4Colour::Yellow(),
                                                   "ds" );

           if (verbosityLevel > 0) {
             G4cout << __func__ << " parent params: " << myPars << G4endl;
             G4cout << __func__ << " core   name:   " << ccrTmpFCore.name << G4endl;
             G4cout << __func__ << " core   params: " << cableRunCoreParams << G4endl;
             // checkForOverlaps( ccrTmpFCore.physical, _config, verbosityLevel>0);
           }

         }

       } // end of if ( CableRunVersion > 1 )

       //Add cabling outside IFB and add panels representing cables exiting IFB Seal
       if ( ds->cableRunVersion() > 2 ) {
         //Defining the outter circular arc for the calorimeter cabling
         TubsParams  calIFBCableRun1Params ( ds->calR1CableRunIFB(),
                                             ds->calR2CableRunIFB(),
                                             ds->zHLCableRunIFB(),
                                             ds->calPhi0CableRunIFB()*CLHEP::degree,
                                             ds->calDPhiCableRunIFB()*CLHEP::degree);

         CLHEP::Hep3Vector calIFBCableRunLoc( 0.0, 0.0, ds->zCCableRunIFB() );

         VolumeInfo icrTmp1 = nestTubs( "CalIFBCableRun1",
                                        calIFBCableRun1Params,
                                        findMaterialOrThrow(ds->materialCalCableRunIFB()),
                                        0,
                                        calIFBCableRunLoc+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                        // dsShieldParent,
                                        parent, // hall air since outside detector
                                        0,
                                        G4Color::Magenta(),
                                        "ds"
                                        );
         // "Fibre Core"
         placeTubeCore ( "CalIFBCableRunCore1",
                         ds->rCableRunCalCoreFract(),
                         ds->rdCableRunCalCoreFract(),
                         ds->dPhiCableRunCalCoreFract(),
                         ds->materialCableRunCalCore(),
                         G4Color::Yellow(),
                         icrTmp1,
                         calIFBCableRun1Params,
                         // "ds",
                         "ds",
                         _config,
                         1
                         );


         TubsParams  calIFBCableRun2Params ( ds->calR1CableRunIFB(),
                                             ds->calR2CableRunIFB(),
                                             ds->zHLCableRunIFB(),
                                             (180.0 - ds->calPhi0CableRunIFB()
                                              - ds->calDPhiCableRunIFB())*CLHEP::degree,
                                             ds->calDPhiCableRunIFB()*CLHEP::degree);

         VolumeInfo icrTmp2 = nestTubs( "CalIFBCableRun2",
                                        calIFBCableRun2Params,
                                        findMaterialOrThrow(ds->materialCalCableRunIFB()),
                                        0,
                                        calIFBCableRunLoc+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                        // dsShieldParent,
                                        parent,
                                        0,
                                        G4Color::Magenta(),
                                        "ds"
                                        );

         // "Fibre Core"
         placeTubeCore ( "CalIFBCableRunCore2",
                         ds->rCableRunCalCoreFract(),
                         ds->rdCableRunCalCoreFract(),
                         ds->dPhiCableRunCalCoreFract(),
                         ds->materialCableRunCalCore(),
                         G4Color::Yellow(),
                         icrTmp2,
                         calIFBCableRun2Params,
                         // "ds",
                         "ds",
                         _config,
                         1
                         );

         //Add radial components of calorimeter outside cabling
         CLHEP::HepRotation * turn = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
         turn->rotateZ(-ds->calPhiECableRunIFB()*CLHEP::degree);
         //add a buffer so box doesn't go beyond the radius of the track cable ring
         double rend0  = ds->calR2CableRunIFB();
         double boxwidth = ds->calEndWCableRunIFB();
         double deltar = sqrt(rend0*rend0 - boxwidth*boxwidth/4.);
         double calIFBCableRunEnd[] = { (deltar - ds->calREndCableRunIFB())/2.,
                                        boxwidth/2.,
                                        ds->zHLCableRunIFB()};
         CLHEP::Hep3Vector calIFBCableRunEndLoc1( (deltar + ds->calREndCableRunIFB())/2.,
                                                 0.0, ds->zCCableRunIFB() );

         calIFBCableRunEndLoc1.rotateZ(ds->calPhiECableRunIFB()*CLHEP::degree);

         VolumeInfo iceTmp1 = nestBox( "CalIFBCableRunEnd1",
                                       calIFBCableRunEnd,
                                       findMaterialOrThrow(ds->materialCalCableRunIFB()),
                                       turn,
                                       calIFBCableRunEndLoc1+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                       // dsShieldParent,
                                       parent,
                                       0,
                                       G4Color::Magenta(),
                                       // "ds"
                                       "ds"
                                       );
         double calIFBCableRunEndCore[] = { calIFBCableRunEnd[0],
                                            calIFBCableRunEnd[1]*ds->rdCableRunCalCoreFract(),
                                            calIFBCableRunEnd[2]*ds->dPhiCableRunCalCoreFract()};

         VolumeInfo icecTmp1 = nestBox( "CalIFBCableRunEndCore1",
                                        calIFBCableRunEndCore,
                                        findMaterialOrThrow(ds->materialCableRunCalCore()),
                                        nullptr,
                                        CLHEP::Hep3Vector(),
                                        iceTmp1,
                                        0,
                                        G4Color::Yellow(),
                                        // "ds"
                                        "ds"
                                        );

         CLHEP::HepRotation * turn2 = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
         turn2->rotateZ((ds->calPhiECableRunIFB()-180.)*CLHEP::degree);

         CLHEP::Hep3Vector calIFBCableRunEndLoc2( (deltar + ds->calREndCableRunIFB())/2.,
                                                 0.0, ds->zCCableRunIFB() );

         calIFBCableRunEndLoc2.rotateZ((180.-ds->calPhiECableRunIFB())*CLHEP::degree);

         VolumeInfo iceTmp2 = nestBox( "CalIFBCableRunEnd2",
                                       calIFBCableRunEnd,
                                       findMaterialOrThrow(ds->materialCalCableRunIFB()),
                                       turn2,
                                       calIFBCableRunEndLoc2+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                       // dsShieldParent,
                                       parent,
                                       0,
                                       G4Color::Magenta(),
                                       // "ds"
                                       "ds"
                                       );

         VolumeInfo icecTmp2 = nestBox( "CalIFBCableRunEndCore2",
                                        calIFBCableRunEndCore,
                                        findMaterialOrThrow(ds->materialCableRunCalCore()),
                                        nullptr,
                                        CLHEP::Hep3Vector(),
                                        iceTmp2,
                                        0,
                                        G4Color::Yellow(),
                                        // "ds"
                                        "ds"
                                        );

         //Cabling at the bottom of the IFB cabling
         double r_bot_ifb     = ds->trkR2CableRunIFB() + 1.; //use bottom of tracker cabling + 1 mm

         // double boxwidth = ds->calEndWCableRunIFB();
         double calIFBCableRunBot[] = { boxwidth/2., //use same box dimensions
                                        ds->calBLCableRunIFB()/2., //length of the piece
                                        ds->zHLCableRunIFB()}; //same width in z as the rest
         CLHEP::Hep3Vector calIFBCableRunBotLoc1( ds->calBCXCableRunIFB(),
                                                 -r_bot_ifb-calIFBCableRunBot[1], ds->zCCableRunIFB() );


         VolumeInfo icbTmp1 = nestBox( "CalIFBCableRunBottom1",
                                       calIFBCableRunBot,
                                       findMaterialOrThrow(ds->materialCalCableRunIFB()),
                                       nullptr,
                                       calIFBCableRunBotLoc1+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       G4Color::Magenta(),
                                       "ds"
                                       );

         double calIFBCableRunBotCore[] = { calIFBCableRunBot[0]*ds->dPhiCableRunCalCoreFract(),
                                            calIFBCableRunBot[1],
                                            calIFBCableRunBot[2]*ds->rdCableRunCalCoreFract()};

         VolumeInfo icbcTmp1 = nestBox( "CalIFBCableRunBottomCore1",
                                        calIFBCableRunBotCore,
                                        findMaterialOrThrow(ds->materialCableRunCalCore()),
                                        nullptr,
                                        CLHEP::Hep3Vector(),
                                        icbTmp1,
                                        0,
                                        G4Color::Yellow(),
                                        "ds"
                                        );

         CLHEP::Hep3Vector calIFBCableRunBotLoc2(-ds->calBCXCableRunIFB(),
                                                 -r_bot_ifb-calIFBCableRunBot[1], ds->zCCableRunIFB() );


         VolumeInfo icbTmp2 = nestBox( "CalIFBCableRunBottom2",
                                       calIFBCableRunBot,
                                       findMaterialOrThrow(ds->materialCalCableRunIFB()),
                                       nullptr,
                                       calIFBCableRunBotLoc2+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       G4Color::Magenta(),
                                       "ds"
                                       );

         VolumeInfo icbcTmp2 = nestBox( "CalIFBCableRunBottomCore2",
                                        calIFBCableRunBotCore,
                                        findMaterialOrThrow(ds->materialCableRunCalCore()),
                                        nullptr,
                                        CLHEP::Hep3Vector(),
                                        icbTmp2,
                                        0,
                                        G4Color::Yellow(),
                                        "ds"
                                        );



         //Define panels on either side of exit of DS representing cables leaving
         TubsParams  calIFBCableExit1Params ( ds->calPR1CableRunIFB(),
                                              ds->calPR2CableRunIFB(),
                                              ds->calPZHLCableRunIFB(),
                                              ds->calPPhi0CableRunIFB()*CLHEP::degree,
                                              ds->calPDPhiCableRunIFB()*CLHEP::degree);

         TubsParams  calIFBCableExit2Params ( ds->calPR1CableRunIFB(),
                                              ds->calPR2CableRunIFB(),
                                              ds->calPZHLCableRunIFB(),
                                              (180.-ds->calPPhi0CableRunIFB()
                                               -ds->calPDPhiCableRunIFB())*CLHEP::degree,
                                              ds->calPDPhiCableRunIFB()*CLHEP::degree);

         CLHEP::Hep3Vector calIFBCablePInLoc( 0.0, 0.0, ds->calPZInCableRunIFB() );
         CLHEP::Hep3Vector calIFBCablePOutLoc( 0.0, 0.0, ds->calPZOutCableRunIFB() );

         VolumeInfo icpTmp1 = nestTubs( "CalIFBCablePanelIn1",
                                        calIFBCableExit1Params,
                                        findMaterialOrThrow(ds->calPMatCableRunIFB()),
                                        0,
                                        calIFBCablePInLoc,
                                        dsShieldPointer,
                                        0,
                                        G4Color::Magenta(),
                                        "ds"
                                        );
         // "Fibre Core"
         placeTubeCore ( "CalIFBCablePanelInCore1",
                         ds->rCableRunCalCoreFract(),
                         ds->rdCableRunCalCoreFract(),
                         ds->dPhiCableRunCalCoreFract(),
                         ds->materialCableRunCalCore(),
                         G4Color::Yellow(),
                         icpTmp1,
                         calIFBCableExit1Params,
                         "ds",
                         _config,
                         0
                         );

         VolumeInfo icpTmp2 = nestTubs( "CalIFBCablePanelPOut1",
                                        calIFBCableExit1Params,
                                        findMaterialOrThrow(ds->calPMatCableRunIFB()),
                                        0,
                                        calIFBCablePOutLoc+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        G4Color::Magenta(),
                                        "ds"
                                        );
         // "Fibre Core"
         placeTubeCore ( "CalIFBCablePanelOutCore1",
                         ds->rCableRunCalCoreFract(),
                         ds->rdCableRunCalCoreFract(),
                         ds->dPhiCableRunCalCoreFract(),
                         ds->materialCableRunCalCore(),
                         G4Color::Yellow(),
                         icpTmp2,
                         calIFBCableExit1Params,
                         "ds",
                         _config,
                         0
                         );
         VolumeInfo icpTmp3 = nestTubs( "CalIFBCablePanelIn2",
                                        calIFBCableExit2Params,
                                        findMaterialOrThrow(ds->calPMatCableRunIFB()),
                                        0,
                                        calIFBCablePInLoc,
                                        dsShieldPointer,
                                        0,
                                        G4Color::Magenta(),
                                        "ds"
                                        );
         // "Fibre Core"
         placeTubeCore ( "CalIFBCablePanelInCore2",
                         ds->rCableRunCalCoreFract(),
                         ds->rdCableRunCalCoreFract(),
                         ds->dPhiCableRunCalCoreFract(),
                         ds->materialCableRunCalCore(),
                         G4Color::Yellow(),
                         icpTmp3,
                         calIFBCableExit2Params,
                         "ds",
                         _config,
                         0
                         );
         VolumeInfo icpTmp4 = nestTubs( "CalIFBCablePanelPOut2",
                                        calIFBCableExit2Params,
                                        findMaterialOrThrow(ds->calPMatCableRunIFB()),
                                        0,
                                        calIFBCablePOutLoc+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        G4Color::Magenta(),
                                        "ds"
                                        );
         // "Fibre Core"
         placeTubeCore ( "CalIFBCablePanelOutCore2",
                         ds->rCableRunCalCoreFract(),
                         ds->rdCableRunCalCoreFract(),
                         ds->dPhiCableRunCalCoreFract(),
                         ds->materialCableRunCalCore(),
                         G4Color::Yellow(),
                         icpTmp4,
                         calIFBCableExit2Params,
                         "ds",
                         _config,
                         0
                         );

         //Define circular arcs for Tracker cabling exiting DS
         TubsParams  trkIFBCableRun1Params ( ds->trkR1CableRunIFB(),
                                             ds->trkR2CableRunIFB(),
                                             ds->zHLCableRunIFB(),
                                             ds->trkPhi0CableRunIFB()*CLHEP::degree,
                                             ds->trkDPhiCableRunIFB()*CLHEP::degree);


         VolumeInfo icrTmp3 = nestTubs( "TrkIFBCableRun1",
                                        trkIFBCableRun1Params,
                                        findMaterialOrThrow(ds->materialTrkCableRunIFB()),
                                        0,
                                        calIFBCableRunLoc+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        G4Color::Magenta(),
                                        "ds"
                                        );

         // "Fibre Core"
         placeTubeCore ( "TrkIFBCableRunCore1",
                         ds->rCableRunTrkCoreFract(),
                         ds->rdCableRunTrkCoreFract(),
                         ds->dPhiCableRunTrkCoreFract(),
                         ds->materialCableRunTrkCore(),
                         G4Color::Yellow(),
                         icrTmp3,
                         trkIFBCableRun1Params,
                         "ds",
                         _config,
                         1
                         );

         TubsParams  trkIFBCableRun2Params ( ds->trkR1CableRunIFB(),
                                             ds->trkR2CableRunIFB(),
                                             ds->zHLCableRunIFB(),
                                             (180.0 - ds->trkPhi0CableRunIFB()
                                              - ds->trkDPhiCableRunIFB())*CLHEP::degree,
                                             ds->trkDPhiCableRunIFB()*CLHEP::degree);

         VolumeInfo icrTmp4 = nestTubs( "TrkIFBCableRun2",
                                        trkIFBCableRun2Params,
                                        findMaterialOrThrow(ds->materialTrkCableRunIFB()),
                                        0,
                                        calIFBCableRunLoc+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        G4Color::Magenta(),
                                        "ds"
                                        );

         // "Fibre Core"
         placeTubeCore ( "TrkIFBCableRunCore2",
                         ds->rCableRunTrkCoreFract(),
                         ds->rdCableRunTrkCoreFract(),
                         ds->dPhiCableRunTrkCoreFract(),
                         ds->materialCableRunTrkCore(),
                         G4Color::Yellow(),
                         icrTmp4,
                         trkIFBCableRun2Params,
                         "ds",
                         _config,
                         1
                         );

         //Define radial components of Tracker cabling exiting DS
         CLHEP::HepRotation * turn3 = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
         turn3->rotateZ(-ds->trkPhiECableRunIFB()*CLHEP::degree);

         double trkIFBCableRunEnd[] = { 0.,
                                        ds->trkEndWCableRunIFB()/2.,
                                        ds->zHLCableRunIFB()};
         //use width and radius of cable runs to find x length to not intesect
         trkIFBCableRunEnd[0] = sqrt(pow(ds->calR1CableRunIFB(),2)-pow(trkIFBCableRunEnd[1],2)) - 5.;
         //half length is half the distance between initial R and intersect point
         trkIFBCableRunEnd[0] =  (trkIFBCableRunEnd[0] - ds->trkREndCableRunIFB())/2.;
         if(trkIFBCableRunEnd[0] < 0.)   trkIFBCableRunEnd[0] = 0.;

         CLHEP::Hep3Vector trkIFBCableRunEndLoc1( ds->trkREndCableRunIFB() + trkIFBCableRunEnd[0],
                                                 0.0, ds->zCCableRunIFB() );
         trkIFBCableRunEndLoc1.rotateZ(ds->trkPhiECableRunIFB()*CLHEP::degree);

         VolumeInfo iceTmp3 = nestBox( "TrkIFBCableRunEnd1",
                                       trkIFBCableRunEnd,
                                       findMaterialOrThrow(ds->materialTrkCableRunIFB()),
                                       turn3,
                                       trkIFBCableRunEndLoc1+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       G4Color::Magenta(),
                                       "ds"
                                       );
         double trkIFBCableRunEndCore[] = { trkIFBCableRunEnd[0],
                                            trkIFBCableRunEnd[1]*ds->rdCableRunTrkCoreFract(),
                                            trkIFBCableRunEnd[2]*ds->dPhiCableRunTrkCoreFract()};

         VolumeInfo icecTmp3 = nestBox( "TrkIFBCableRunEndCore1",
                                        trkIFBCableRunEndCore,
                                        findMaterialOrThrow(ds->materialCableRunTrkCore()),
                                        nullptr,
                                        CLHEP::Hep3Vector(),
                                        iceTmp3,
                                        0,
                                        G4Color::Yellow(),
                                        "ds"
                                        );

         CLHEP::HepRotation * turn4 = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
         turn4->rotateZ((ds->trkPhiECableRunIFB()-180.)*CLHEP::degree);

         CLHEP::Hep3Vector trkIFBCableRunEndLoc2( ds->trkREndCableRunIFB() + trkIFBCableRunEnd[0],
                                                 0.0, ds->zCCableRunIFB() );
         trkIFBCableRunEndLoc2.rotateZ((180.-ds->trkPhiECableRunIFB())*CLHEP::degree);

         VolumeInfo iceTmp4 = nestBox( "TrkIFBCableRunEnd2",
                                       trkIFBCableRunEnd,
                                       findMaterialOrThrow(ds->materialTrkCableRunIFB()),
                                       turn4,
                                       trkIFBCableRunEndLoc2+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       G4Color::Magenta(),
                                       "ds"
                                       );
         VolumeInfo icecTmp4 = nestBox( "TrkIFBCableRunEndCore2",
                                        trkIFBCableRunEndCore,
                                        findMaterialOrThrow(ds->materialCableRunTrkCore()),
                                        nullptr,
                                        CLHEP::Hep3Vector(),
                                        iceTmp4,
                                        0,
                                        G4Color::Yellow(),
                                        "ds"
                                        );


         //Cabling at the bottom of the IFB cabling
         double trkIFBCableRunBot[] = { ds->trkEndWCableRunIFB()/2., //use same box dimensions
                                        ds->trkBLCableRunIFB()/2., //length of the piece
                                        ds->zHLCableRunIFB()}; //same width in z as the rest
         CLHEP::Hep3Vector trkIFBCableRunBotLoc1( ds->trkBCXCableRunIFB(),
                                                 -r_bot_ifb-trkIFBCableRunBot[1], ds->zCCableRunIFB() );
         VolumeInfo icbTmp3 = nestBox( "TrkIFBCableRunBottom1",
                                       trkIFBCableRunBot,
                                       findMaterialOrThrow(ds->materialTrkCableRunIFB()),
                                       nullptr,
                                       trkIFBCableRunBotLoc1+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       G4Color::Magenta(),
                                       "ds"
                                       );

         double trkIFBCableRunBotCore[] = { trkIFBCableRunBot[0]*ds->dPhiCableRunTrkCoreFract(),
                                            trkIFBCableRunBot[1],
                                            trkIFBCableRunBot[2]*ds->rdCableRunTrkCoreFract()};

         VolumeInfo icbcTmp3 = nestBox( "TrkIFBCableRunBottomCore1",
                                        trkIFBCableRunBotCore,
                                        findMaterialOrThrow(ds->materialCableRunTrkCore()),
                                        nullptr,
                                        CLHEP::Hep3Vector(),
                                        icbTmp3,
                                        0,
                                        G4Color::Yellow(),
                                        "ds"
                                        );

         CLHEP::Hep3Vector trkIFBCableRunBotLoc2(-ds->trkBCXCableRunIFB(),
                                                 -r_bot_ifb-trkIFBCableRunBot[1], ds->zCCableRunIFB() );

         VolumeInfo icbTmp4 = nestBox( "TrkIFBCableRunBottom2",
                                       trkIFBCableRunBot,
                                       findMaterialOrThrow(ds->materialTrkCableRunIFB()),
                                       nullptr,
                                       trkIFBCableRunBotLoc2+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       G4Color::Magenta(),
                                       "ds"
                                       );

         VolumeInfo icbcTmp4 = nestBox( "TrkIFBCableRunBottomCore2",
                                        trkIFBCableRunBotCore,
                                        findMaterialOrThrow(ds->materialCableRunTrkCore()),
                                        nullptr,
                                        CLHEP::Hep3Vector(),
                                        icbTmp4,
                                        0,
                                        G4Color::Yellow(),
                                        "ds"
                                        );


         //Define panels on either side of exit of DS representing cables leaving for tracker
         TubsParams  trkIFBCableExit1Params ( ds->trkPR1CableRunIFB(),
                                              ds->trkPR2CableRunIFB(),
                                              ds->trkPZHLCableRunIFB(),
                                              ds->trkPPhi0CableRunIFB()*CLHEP::degree,
                                              ds->trkPDPhiCableRunIFB()*CLHEP::degree);
         TubsParams  trkIFBCableExit2Params ( ds->trkPR1CableRunIFB(),
                                              ds->trkPR2CableRunIFB(),
                                              ds->trkPZHLCableRunIFB(),
                                              (180.-ds->trkPPhi0CableRunIFB()
                                               -ds->trkPDPhiCableRunIFB())*CLHEP::degree,
                                              ds->trkPDPhiCableRunIFB()*CLHEP::degree);

         CLHEP::Hep3Vector trkIFBCablePInLoc( 0.0, 0.0, ds->trkPZInCableRunIFB() );
         CLHEP::Hep3Vector trkIFBCablePOutLoc( 0.0, 0.0, ds->trkPZOutCableRunIFB() );

         VolumeInfo icpTmp5 = nestTubs( "TrkIFBCablePanelIn1",
                                        trkIFBCableExit1Params,
                                        findMaterialOrThrow(ds->trkPMatCableRunIFB()),
                                        0,
                                        trkIFBCablePInLoc,
                                        dsShieldPointer,
                                        0,
                                        G4Color::Magenta(),
                                        "ds"
                                        );

         // "Fibre Core"
         placeTubeCore ( "TrkIFBCablePanelInCore1",
                         ds->rCableRunTrkCoreFract(),
                         ds->rdCableRunTrkCoreFract(),
                         ds->dPhiCableRunTrkCoreFract(),
                         ds->materialCableRunTrkCore(),
                         G4Color::Yellow(),
                         icpTmp5,
                         trkIFBCableExit1Params,
                         "ds",
                         _config,
                         0
                         );

         VolumeInfo icpTmp6 = nestTubs( "TrkIFBCablePanelPOut1",
                                        trkIFBCableExit1Params,
                                        findMaterialOrThrow(ds->trkPMatCableRunIFB()),
                                        0,
                                        trkIFBCablePOutLoc+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        G4Color::Magenta(),
                                        "ds"
                                        );

         // "Fibre Core"
         placeTubeCore ( "TrkIFBCablePanelOutCore1",
                         ds->rCableRunTrkCoreFract(),
                         ds->rdCableRunTrkCoreFract(),
                         ds->dPhiCableRunTrkCoreFract(),
                         ds->materialCableRunTrkCore(),
                         G4Color::Yellow(),
                         icpTmp6,
                         trkIFBCableExit1Params,
                         "ds",
                         _config,
                         0
                         );

         VolumeInfo icpTmp7 = nestTubs( "TrkIFBCablePanelIn2",
                                        trkIFBCableExit2Params,
                                        findMaterialOrThrow(ds->trkPMatCableRunIFB()),
                                        0,
                                        trkIFBCablePInLoc,
                                        dsShieldPointer,
                                        0,
                                        G4Color::Magenta(),
                                        "ds"
                                        );

         // "Fibre Core"
         placeTubeCore ( "TrkIFBCablePanelInCore2",
                         ds->rCableRunTrkCoreFract(),
                         ds->rdCableRunTrkCoreFract(),
                         ds->dPhiCableRunTrkCoreFract(),
                         ds->materialCableRunTrkCore(),
                         G4Color::Yellow(),
                         icpTmp7,
                         trkIFBCableExit2Params,
                         "ds",
                         _config,
                         0
                         );

         VolumeInfo icpTmp8 = nestTubs( "TrkIFBCablePanelPOut2",
                                        trkIFBCableExit2Params,
                                        findMaterialOrThrow(ds->trkPMatCableRunIFB()),
                                        0,
                                        trkIFBCablePOutLoc+dsShieldParent.centerInMu2e()-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        G4Color::Magenta(),
                                        "ds"
                                        );

         // "Fibre Core"
         placeTubeCore ( "TrkIFBCablePanelOutCore2",
                         ds->rCableRunTrkCoreFract(),
                         ds->rdCableRunTrkCoreFract(),
                         ds->dPhiCableRunTrkCoreFract(),
                         ds->materialCableRunTrkCore(),
                         G4Color::Yellow(),
                         icpTmp8,
                         trkIFBCableExit2Params,
                         "ds",
                         _config,
                         0
                         );


       } // end of if ( CableRunVersion > 2 )
     } // end of if ( ds->hasCableRunCal() )

     if ( ds->hasCableRunTrk() ) {

       TubsParams  trkCableRun1Params ( ds->rInCableRunTrk(),
                                        ds->rOutCableRunTrk(),
                                        ds->lengthCableRunTrk(),
                                        ds->phi0CableRunTrk()*CLHEP::degree,
                                        ds->dPhiCableRunTrk()*CLHEP::degree);

       CLHEP::Hep3Vector trkCableRunLoc( 0.0, 0.0, ds->zCCableRunTrk() );

       VolumeInfo tcrTmp1 = nestTubs( "TrkCableRun1",
                                      trkCableRun1Params,
                                      findMaterialOrThrow(ds->trkCableRunMaterial()),
                                      0,
                                      trkCableRunLoc,
                                      dsShieldPointer,
                                      0,
                                      G4Color::Magenta(),
                                      "ds"
                                      );

       if ( ds->cableRunVersion() > 2 ) {

         // "Fibre Core"
         placeTubeCore ( "TrkCableRun1Core",
                         ds->rCableRunTrkCoreFract(),
                         ds->rdCableRunTrkCoreFract(),
                         ds->dPhiCableRunCalCoreFract(),
                         ds->materialCableRunTrkCore(),
                         G4Color::Yellow(),
                         tcrTmp1,
                         trkCableRun1Params,
                         "ds",
                         _config
                         );

       }

       // Now the second one
       TubsParams  trkCableRun2Params ( ds->rInCableRunTrk(),
                                        ds->rOutCableRunTrk(),
                                        ds->lengthCableRunTrk(),
                                        (180.0 - ds->phi0CableRunTrk()
                                         - ds->dPhiCableRunTrk())
                                        *CLHEP::degree,
                                        ds->dPhiCableRunTrk()*CLHEP::degree);

       VolumeInfo tcrTmp2=nestTubs( "TrkCableRun2",
                                    trkCableRun2Params,
                                    findMaterialOrThrow(ds->trkCableRunMaterial()),
                                    0,
                                    trkCableRunLoc,
                                    dsShieldPointer,
                                    0,
                                    G4Color::Magenta(),
                                    "ds"
                                    );

       if ( ds->cableRunVersion() > 2 ) {

         // "Fibre Core"
         placeTubeCore ( "TrkCableRun2Core",
                         ds->rCableRunTrkCoreFract(),
                         ds->rdCableRunTrkCoreFract(),
                         ds->dPhiCableRunCalCoreFract(),
                         ds->materialCableRunTrkCore(),
                         G4Color::Yellow(),
                         tcrTmp2,
                         trkCableRun2Params,
                         "ds",
                         _config
                         );

       }

     } // end of if ( ds->hasCableRunTrk() )

     if ( ds->hasServicePipes() ) {
       TubsParams pipeMomParms( 0.0,
                                ds->servicePipeROut(),
                                ds->servicePipeHalfLength());//default phi,dphi
       TubsParams pipeSelfParms(ds->servicePipeRIn(),
                                ds->servicePipeROut(),
                                ds->servicePipeHalfLength() );
       std::vector<double>theXs = ds->servicePipeXCs();
       double theY = ds->servicePipeYC();
       double theZ = ds->servicePipeZC();

       for ( unsigned int iP = 0; iP < theXs.size(); iP++ ) {
         CLHEP::Hep3Vector pipeLoc(theXs[iP],theY, theZ);
         ostringstream pmName;
         pmName << "DSservicePipeMother" << iP+1;
         VolumeInfo DSPipeMother = nestTubs ( pmName.str(),
                                              pipeMomParms,
                                              findMaterialOrThrow(ds->servicePipeFillMat()),
                                              0,
                                              pipeLoc,
                                              dsShieldParent,
                                              0,
                                              G4Color::Magenta(),
                                              "ds"
                                              );

         // Now make it a real pipe
         ostringstream pName;
         pName << "DSservicePipe" << iP+1;
         nestTubs ( pName.str(),
                    pipeSelfParms,
                    findMaterialOrThrow(ds->servicePipeMaterial()),
                    0,
                    CLHEP::Hep3Vector(),
                    DSPipeMother,
                    0,
                    G4Color::Magenta(),
                    "ds"
                    );

       } // End of loop over service pipes
     } // end of if ( ds->hasServicePipes() )

  } // end of Mu2eWorld::constructDS;

  // limited utility function to place the core tube in a tube
  void placeTubeCore ( const std::string & name,
                       double radiusFract,
                       double radiusDFract,
                       double dPhiFraction,
                       const std::string & material,
                       const G4Colour & color,
                       const VolumeInfo & parent,
                       const TubsParams & parentParams,
                       const std::string &  lookupToken,
                       const SimpleConfig & config,
                       const int zNotPhi
                       ) {

    int const verbosityLevel = config.getInt("ds.verbosityLevel",0);
    TubsParams cableRunCoreParams =
      calculateTubeCoreParams(parentParams,
                              radiusFract,
                              radiusDFract,
                              dPhiFraction,
                              verbosityLevel,
                              zNotPhi);

    VolumeInfo tempCore = nestTubs( name,
                                    cableRunCoreParams,
                                    findMaterialOrThrow(material),
                                    nullptr,
                                    CLHEP::Hep3Vector(),
                                    parent,
                                    0,
                                    color,
                                    lookupToken
                                    );
    if (verbosityLevel > 0) {
      G4cout << __func__ << " parent params: " << parentParams << G4endl;
      G4cout << __func__ << " core   name:   " << name << G4endl;
      G4cout << __func__ << " core   params: " << cableRunCoreParams << G4endl;
      // checkForOverlaps( tempCore.physical, config, verbosityLevel>0);
    }

  }

  TubsParams calculateTubeCoreParams (const TubsParams& parentParams,
                                      double radiusFract,
                                      double radiusDFract,
                                      double dPhiFraction,
                                      int verbosityLevel,
                                      const int zNotPhi) {


    // calculating the core parameters based on the the parent tube
    // sizes and the core parameters relative to the tube itself

    // radial

    double coreCenterRadius = parentParams.outerRadius() -
      ( parentParams.outerRadius() - parentParams.innerRadius() ) * radiusFract;
    // 0 -> rOut; 1 -> rIn; 0.5 -> (rIn + rOut)/2;

    double coreRadialHalfExtent =
      ( parentParams.outerRadius() - parentParams.innerRadius() ) * 0.5 *
      radiusDFract;

    double coreInnerRadius = coreCenterRadius - coreRadialHalfExtent;
    double coreOuterRadius = coreCenterRadius + coreRadialHalfExtent;

    if ( coreInnerRadius < parentParams.innerRadius() ||
         coreOuterRadius > parentParams.outerRadius() ) {

      throw cet::exception("GEOM") << __func__
                                   << " inconsitent cable core parameters: "
                                   << coreInnerRadius
                                   << ", "
                                   << coreOuterRadius
                                   << "\n";

    }

    // angular (parent angular params are multiplied by CLHEP::degree already)

    if(zNotPhi == 0) { //if a small arc that extends in z plane, with dPhi and dR fractions
      double coreCenterPhi = parentParams.phi0() + 0.5 * parentParams.phiTotal();
      double coreAngularHalfExtent = parentParams.phiTotal() * dPhiFraction * 0.5;
      double corePhi0 = coreCenterPhi - coreAngularHalfExtent;
      double coreDeltaPhi = 2.0 * coreAngularHalfExtent;

      if ( corePhi0 < parentParams.phi0() ) {

        throw cet::exception("GEOM") << __func__
                                     << " inconsitent cable core parameters: "
                                     << corePhi0
                                     << ", "
                                     << coreDeltaPhi
                                     << "\n";

      }
      return TubsParams( coreInnerRadius,
                         coreOuterRadius,
                         parentParams.zHalfLength(),
                         corePhi0,
                         coreDeltaPhi );

    } else { //if an arc in a z plane, with dZ and dR fractions
      double coreZHalfExtent = parentParams.zHalfLength() * dPhiFraction;

      if ( dPhiFraction > 1. || dPhiFraction < 0. ) {

        throw cet::exception("GEOM") << __func__
                                     << " inconsitent cable core parameters: "
                                     << coreZHalfExtent
                                     << ", "
                                     << dPhiFraction
                                     << "\n";

      }
      return TubsParams( coreInnerRadius,
                         coreOuterRadius,
                         coreZHalfExtent,
                         parentParams.phi0(),
                         parentParams.phiTotal());

    }
  }


}
