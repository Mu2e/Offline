//
// Free function to create Transport Solenoid
//
//
// Original author KLG based on Mu2eWorld constructTS
//
// Notes:
// Construct the TS.  Parent volume is the air inside of the hall.

// Mu2e includes.
#include "Offline/Mu2eG4/inc/constructTS.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"

// C++ includes
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes.
#include "Offline/BeamlineGeom/inc/Collimator_TS1.hh"
#include "Offline/BeamlineGeom/inc/PbarWindow.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/GeomPrimitives/inc/PolyconsParams.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/nestPolycone.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestTorus.hh"
#include "Offline/Mu2eG4/inc/nestCons.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/GeomPrimitives/inc/Polycone.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSShield.hh"

// G4 includes
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Cons.hh"
#include "Geant4/G4Torus.hh"
#include "Geant4/G4Trd.hh"
#include "Geant4/G4IntersectionSolid.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4TwoVector.hh"
#include "Geant4/G4ExtrudedSolid.hh"
#include "Geant4/G4GeometryManager.hh"


using namespace std;

namespace mu2e {

  void constructTS( VolumeInfo const & parent,
                    SimpleConfig const & c
                    ) {

    GeomHandle<Beamline> bl;
    Mu2eG4Helper* _helper = &(*art::ServiceHandle<Mu2eG4Helper>() );

    constructCryostat   ( parent, c, *bl);
    constructCoils      ( parent, c, *bl);
    constructCAs        ( parent, c, *bl);
    constructCollimators( parent, c, *bl);
    constructDegrader   ( _helper->locateVolInfo("TS5Vacuum"), c, *bl);
    constructPbarWindow ( _helper->locateVolInfo("TS3Vacuum"), c, *bl);

  }

  //__________________________________
  //
  // CONSTRUCT CRYOSTAT
  //__________________________________

  void constructCryostat( VolumeInfo const& parent,
                          SimpleConfig const& config,
                          Beamline const& bl ) {

    TransportSolenoid const * ts     ( &bl.getTS() );
    StraightSection   const * strsec (nullptr);
    TorusSection      const * torsec (nullptr);

    GeomHandle<PSShield> hrs;

    Mu2eG4Helper * const _helper = &(*art::ServiceHandle<Mu2eG4Helper>());

    const int  verbosityLevel      = config.getInt ("ts.cryo.verbosityLevel", 0     );
    const bool polyLiningUp        = config.getBool("ts.polyliner.Up.build"   , false );
    const bool polyLiningDown      = config.getBool("ts.polyliner.Down.build" , false );
    const bool tsThermalShield     = config.getBool("ts.thermalshield.build"  , false );

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "TSCryo", "ts.cryo"      );
    geomOptions->loadEntry( config, "TSPoly", "ts.polyliner" );
    geomOptions->loadEntry( config, "TSCA"  , "ts.cas"       );
    geomOptions->loadEntry( config, "PSShield", "PSShield" );
    // For how all pieces are made from one of two types of material,
    // vacuum or average coils + cryostat material.
    G4Material* downstreamVacuumMaterial  = findMaterialOrThrow(ts->downstreamVacuumMaterial());
    G4Material* upstreamVacuumMaterial    = findMaterialOrThrow(ts->upstreamVacuumMaterial());
    G4Material* cryoMaterial              = findMaterialOrThrow(ts->material());
    G4Material* thermalShieldMLIMaterial  = (tsThermalShield) ? findMaterialOrThrow(ts->thermalShieldMLIMaterial()) : 0;
    G4Material* thermalShieldMidMaterial  = (tsThermalShield) ? findMaterialOrThrow(ts->thermalShieldMidMaterial()) : 0;

    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

    // Build upstream end wall of TS1
    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS1,TransportSolenoid::TSRadialPart::IN );

    CLHEP::Hep3Vector pos( strsec->getGlobal().x(),
                           strsec->getGlobal().y(),
                           strsec->getGlobal().z()-strsec->getHalfLength()+ts->endWallU1_halfLength() );

    std::string tssName  = "TS1UpstreamEndwall";

    G4VSolid *tssolid= new G4Tubs( "TS1UpstreamEndwallFull", ts->endWallU1_rIn(),ts->endWallU1_rOut(),ts->endWallU1_halfLength(),0.,2*M_PI);

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    // Get the beam pipe
    G4Tubs* beamPassTub = nullptr;
    CLHEP::HepRotation* turn = nullptr;
    CLHEP::Hep3Vector place;
    getBeamPipe(config, reg, beamPassTub, turn, place);

    G4SubtractionSolid* aSolid =
      new G4SubtractionSolid ( "TS1UpstreamEndwallFull",
                               tssolid,
                               beamPassTub,
                               turn,
                               place-pos);

    VolumeInfo tss(tssName,
                   pos-_hallOriginInMu2e,
                   parent.centerInWorld);
    tss.solid = aSolid;
    tss.solid->SetName(tss.name);

    finishNesting(tss,
                  cryoMaterial,
                  0,
                  tss.centerInParent,
                  parent.logical,
                  0,
                  true,
                  //G4Colour(0xFF/double(0xFF), 0x99/double(0xFF), 0),
                  G4Color::Blue(),
                  true,
                  true,
                  true,
                  false
                  );
    // Build TS1
    strsec = ts->getTSVacuum<StraightSection>(TransportSolenoid::TSRegion::TS1);
    nestTubs( "TS1Vacuum",
              TubsParams( strsec->rIn(),
                          strsec->rOut(),
                          strsec->getHalfLength() ),
              upstreamVacuumMaterial,
              strsec->getRotation(),
              strsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );

    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS1,TransportSolenoid::TSRadialPart::IN);
    tssName = "TS1InnerCryoShell";
    double ts1InsVacRIn = strsec->rOut();
    nestTubs( tssName,
              TubsParams( strsec->rIn(),
                          strsec->rOut(),
                          strsec->getHalfLength() ),
              cryoMaterial,
              strsec->getRotation(),
              strsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );

    verbosityLevel &&
      G4cout << __func__ << " " << tssName << " Mass in kg: "
             << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg
             << G4endl;

    if ( verbosityLevel > 0) {
      G4cout << __func__ << " TS1(in)  OffsetInMu2e  : " << strsec->getGlobal()   << G4endl;
      G4cout << __func__ << " TS1(in)  Extent        :[ " << strsec->getGlobal().z() - strsec->getHalfLength() <<","
             << strsec->getGlobal().z() + strsec->getHalfLength() << "]" << G4endl;
      G4cout << __func__ << " TS1(in)  rotation      : " << *(strsec->getRotation()) << G4endl;
    }

    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS1,TransportSolenoid::TSRadialPart::OUT );
    tssName =  "TS1OuterCryoShell";
    double ts1InsVacROut = strsec->rIn();

    G4VSolid *tscryosolid= new G4Tubs( "TS1OuterCryoShellFull", strsec->rIn(),strsec->rOut(),strsec->getHalfLength(),0.,2*M_PI);

    G4SubtractionSolid* acryoSolid =
      new G4SubtractionSolid ( "TS1OuterCryoShell",
                               tscryosolid,
                               beamPassTub,
                               turn,
                               place-strsec->getGlobal());
    VolumeInfo tssCryo(tssName,
                       pos-_hallOriginInMu2e,
                       parent.centerInWorld);
    tssCryo.solid = acryoSolid;
    tssCryo.solid->SetName(tssCryo.name);

    finishNesting(tssCryo,
                  cryoMaterial,
                  strsec->getRotation(),
                  strsec->getGlobal()-_hallOriginInMu2e,
                  parent.logical,
                  0,
                  true,
                  //G4Colour(0xFF/double(0xFF), 0x99/double(0xFF), 0),
                  G4Color::Blue(),
                  true,
                  true,
                  true,
                  false
                  );

    // Build downstream partial end wall of TS1
    CLHEP::Hep3Vector pos2( strsec->getGlobal().x(),
                            strsec->getGlobal().y(),
                            strsec->getGlobal().z()+strsec->getHalfLength()+ts->endWallU2_halfLength() );

    tssName = "TS1DownstreamEndwall";
    G4VSolid *tsdsolid= new G4Tubs( "TS1DownstreamEndwallFull", ts->endWallU2_rIn(),ts->endWallU2_rOut(),ts->endWallU2_halfLength(),0.,2*M_PI);


    G4SubtractionSolid* adSolid =
      new G4SubtractionSolid ( "TS1DownstreamEndwallFull",
                               tsdsolid,
                               beamPassTub,
                               turn,
                               place-pos2);

    VolumeInfo tsd(tssName,
                   pos2-_hallOriginInMu2e,
                   parent.centerInWorld);
    tsd.solid = adSolid;
    tsd.solid->SetName(tsd.name);

    finishNesting(tsd,
                  cryoMaterial,
                  0,
                  pos2-_hallOriginInMu2e,
                  parent.logical,
                  0,
                  true,
                  G4Color::Blue(),
                  true,
                  true,
                  true,
                  false
                  );

    // Put in the insulating vacuum, which will serve as the mother volume
    // for the coils and coil assemblies (CAs).
    double ts1InsVacHalfLen = strsec->getHalfLength() - ts->endWallU1_halfLength() + ts->endWallU2_halfLength();
    double ts1InsVacOffset = ts->endWallU1_halfLength() + ts->endWallU2_halfLength();
    CLHEP::Hep3Vector ts1CIVLoc = strsec->getGlobal()-_hallOriginInMu2e + CLHEP::Hep3Vector(0,0,ts1InsVacOffset);
    tssName =  "TS1CryoInsVac";
    nestTubs( tssName,
              TubsParams( ts1InsVacRIn,
                          ts1InsVacROut,
                          ts1InsVacHalfLen ),
              upstreamVacuumMaterial,
              strsec->getRotation(),
              ts1CIVLoc,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );

    if(tsThermalShield) {
      VolumeInfo useAsParent = _helper->locateVolInfo( "TS1CryoInsVac" );
      addThermalShield(*ts, useAsParent, config, TransportSolenoid::TSRegion::TS1,
                       thermalShieldMLIMaterial, thermalShieldMidMaterial, -1.);
    }

    // Build TS2
    torsec = ts->getTSVacuum<TorusSection>(TransportSolenoid::TSRegion::TS2);
    nestTorus("TS2Vacuum",
              torsec->getParameters(),
              upstreamVacuumMaterial,
              torsec->getRotation(),
              torsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Green(),
              "TSCryo"
              );

    torsec = ts->getTSPolyLining(TransportSolenoid::TSRegion::TS2);
    if ( polyLiningUp && torsec->rIn() > 0. ) {
      VolumeInfo ts2Vacuum = art::ServiceHandle<Mu2eG4Helper>()->locateVolInfo("TS2Vacuum");
      nestTorus("TS2PolyLining",
                torsec->getParameters(),
                findMaterialOrThrow(torsec->getMaterial()),
                torsec->getRotation(),
                torsec->getGlobal()-ts2Vacuum.centerInMu2e(),
                ts2Vacuum,
                0,
                G4Color::Yellow(),
                "TSPoly"
                );
    }

    torsec = ts->getTSCryo<TorusSection>(TransportSolenoid::TSRegion::TS2,TransportSolenoid::TSRadialPart::IN );
    std::array<double,5> ts2Cryo1Params { { torsec->rIn(), torsec->rOut(), torsec->torusRadius(), torsec->phiStart(), torsec->deltaPhi() } };

    tssName = "TS2InnerCryoShell";
    double ts2InsVacRIn = torsec->rOut();
    nestTorus(tssName,
              ts2Cryo1Params,
              cryoMaterial,
              torsec->getRotation(),
              torsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );

    verbosityLevel &&
      G4cout << __func__ << " " << tssName << " Mass in kg: "
             << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg
             << G4endl;

    torsec = ts->getTSCryo<TorusSection>(TransportSolenoid::TSRegion::TS2,TransportSolenoid::TSRadialPart::OUT );

    tssName = "TS2OuterCryoShell";
    double ts2InsVacROut = torsec->rIn();

    G4VSolid *ts2cryosolid= new G4Torus( "TS2OuterCryoShellFull",torsec->rIn(),torsec->rOut(), torsec->torusRadius(), torsec->phiStart(), torsec->deltaPhi());
    CLHEP::Hep3Vector pos2t=torsec->getGlobal();

    CLHEP::HepRotation const RA= *(torsec->getRotation());
    CLHEP::HepRotation const RAinv = torsec->getRotation()->inverse();
    CLHEP::HepRotation const& RB = *turn;
    CLHEP::HepRotation R_relative_matrix = RB*RAinv;
    CLHEP::HepRotation* R_relative_ptr = reg.add(R_relative_matrix);
    CLHEP::Hep3Vector rotated_place= RA*(place-pos2t);
    G4SubtractionSolid* acryo2Solid =
      new G4SubtractionSolid ( "TS2OuterCryoShell",
                               ts2cryosolid,
                               beamPassTub,
                               R_relative_ptr,
                               rotated_place);
    VolumeInfo ts2Cryo(tssName,
                       pos2t-_hallOriginInMu2e,
                       parent.centerInWorld);
    ts2Cryo.solid = acryo2Solid;
    ts2Cryo.solid->SetName(ts2Cryo.name);

    finishNesting(ts2Cryo,
                  cryoMaterial,
                  torsec->getRotation(),
                  torsec->getGlobal()-_hallOriginInMu2e,
                  parent.logical,
                  0,
                  true,
                  G4Color::Blue(),
                  true,
                  true,
                  true,
                  false
                  );
    /* This is for debugging: create the beamPass volume to see it even if it creates interferences with other mother volumes
       std::string bTubeName  = "beamPass";
       nestTubs( bTubeName,
       TubsParams(  beamPassTub->GetInnerRadius(),
       beamPassTub->GetOuterRadius(),
       beamPassTub->GetZHalfLength() ),
       cryoMaterial,
       turn,
       place- parent.centerInWorld,
       parent,
       0,
       G4Color::Red(),
       "TSCryo"
       );
    */
    // Put in the insulating vacuum, which will serve as the mother volume
    // for the coils and coil assemblies (CAs).
    std::array<double,5> ts2CIVParams { { ts2InsVacRIn, ts2InsVacROut, torsec->torusRadius(), torsec->phiStart(), torsec->deltaPhi() } };


    tssName = "TS2CryoInsVac";
    nestTorus(tssName,
              ts2CIVParams,
              upstreamVacuumMaterial,
              torsec->getRotation(),
              torsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );

    if(tsThermalShield) {
      VolumeInfo useAsParent = _helper->locateVolInfo( "TS2CryoInsVac" );
      addThermalShield(*ts, useAsParent, config, TransportSolenoid::TSRegion::TS2,
                       thermalShieldMLIMaterial, thermalShieldMidMaterial, -1.);
    }

    // Build TS3
    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS3,TransportSolenoid::TSRadialPart::IN );
    nestTubs( "TS3Vacuum",
              TubsParams( 0., ts->ts3InnerRadius(), strsec->getHalfLength() ),
              downstreamVacuumMaterial,
              strsec->getRotation(),
              strsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Green(),
              "TSCryo"
              );

    double ts3InsVacRIn = strsec->rOut();
    tssName = "TS3InnerCryoShell";
    nestTubs( tssName,
              TubsParams( strsec->rIn(),
                          strsec->rOut(),
                          strsec->getHalfLength() ),
              cryoMaterial,
              strsec->getRotation(),
              strsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );

    verbosityLevel &&
      G4cout << __func__ << " " << tssName << " Mass in kg: "
             << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg
             << G4endl;

    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS3,TransportSolenoid::TSRadialPart::OUT );

    tssName = "TS3OuterCryoShell";
    double ts3InsVacROut = strsec->rIn();
    nestTubs( tssName,
              TubsParams( strsec->rIn(),
                          strsec->rOut(),
                          strsec->getHalfLength() ),
              cryoMaterial,
              strsec->getRotation(),
              strsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );

    verbosityLevel &&
      G4cout << __func__ << " " << tssName << " Mass in kg: "
             << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg
             << G4endl;

    if ( verbosityLevel > 0) {
      G4cout << __func__ << " TS3  OffsetInMu2e : " << strsec->getGlobal()   << G4endl;
      G4cout << __func__ << " TS3  rotation     : " << *(strsec->getRotation()) << G4endl;
    }

    // Put in the insulating vacuum, which will serve as the mother volume
    // for the coils and coil assemblies (CAs).

    tssName =  "TS3CryoInsVac";
    nestTubs( tssName,
              TubsParams( ts3InsVacRIn,
                          ts3InsVacROut,
                          strsec->getHalfLength() ),
              upstreamVacuumMaterial,
              strsec->getRotation(),
              strsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );

    if(tsThermalShield) {
      //split each volume into an upstream and downstream end
      double ts3WallThickness = config.getDouble("pbar.support.midThickness")*CLHEP::mm;
      VolumeInfo useAsParent = _helper->locateVolInfo( "TS3CryoInsVac" );
      addThermalShield(*ts, useAsParent, config, TransportSolenoid::TSRegion::TS3,
                       thermalShieldMLIMaterial, thermalShieldMidMaterial, ts3WallThickness);
    }

    //Add TSu-TSd interconnect
    if(config.getBool("ts.tsudinterconnect.build", false)) {
      //strsec = TS3 section
      std::string tsudinterconnect_material = config.getString("ts.tsudinterconnect.material");
      double tsudinterconnect_xoffset = config.getDouble("ts.tsudinterconnect.xoffset");
      double tsudinterconnect_rin = config.getDouble("ts.tsudinterconnect.rIn");
      double tsudinterconnect_rout = config.getDouble("ts.tsudinterconnect.rOut");
      double tsudinterconnect_halflength = config.getDouble("ts.tsudinterconnect.halfLength");
      G4Material* tsud_mat = findMaterialOrThrow(tsudinterconnect_material);

      nestTubs( "TSudInterconnectu",
                TubsParams( tsudinterconnect_rin,
                            tsudinterconnect_rout,
                            tsudinterconnect_halflength),
                tsud_mat,
                strsec->getRotation(),
                strsec->getGlobal()-_hallOriginInMu2e
                + CLHEP::Hep3Vector(tsudinterconnect_xoffset+tsudinterconnect_halflength, 0., 0.),
                parent,
                0,
                G4Color::Red(),
                "TSCryo"
                );
      nestTubs( "TSudInterconnectd",
                TubsParams( tsudinterconnect_rin,
                            tsudinterconnect_rout,
                            tsudinterconnect_halflength),
                tsud_mat,
                strsec->getRotation(),
                strsec->getGlobal()-_hallOriginInMu2e
                - CLHEP::Hep3Vector(tsudinterconnect_xoffset+tsudinterconnect_halflength, 0., 0.),
                parent,
                0,
                G4Color::Red(),
                "TSCryo"
                );
    }

    // Build TS4
    torsec = ts->getTSVacuum<TorusSection>(TransportSolenoid::TSRegion::TS4);
    nestTorus("TS4Vacuum",
              torsec->getParameters(),
              downstreamVacuumMaterial,
              torsec->getRotation(),
              torsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Yellow(),
              "TSCryo"
              );

    torsec = ts->getTSPolyLining(TransportSolenoid::TSRegion::TS4);
    if ( polyLiningDown && torsec->rIn() > 0. ) {
      VolumeInfo ts4Vacuum = art::ServiceHandle<Mu2eG4Helper>()->locateVolInfo("TS4Vacuum");
      nestTorus("TS4PolyLining",
                torsec->getParameters(),
                findMaterialOrThrow( torsec->getMaterial() ),
                torsec->getRotation(),
                torsec->getGlobal()-ts4Vacuum.centerInMu2e(),
                ts4Vacuum,
                0,
                G4Color::Yellow(),
                "TSPoly"
                );
    }

    torsec = ts->getTSCryo<TorusSection>(TransportSolenoid::TSRegion::TS4,TransportSolenoid::TSRadialPart::IN );
    std::array<double,5> ts4Cryo1Params { { torsec->rIn(), torsec->rOut(), torsec->torusRadius(), torsec->phiStart(), torsec->deltaPhi() } };
    tssName = "TS4InnerCryoShell";
    double ts4InsVacRIn = torsec->rOut();
    nestTorus(tssName,
              ts4Cryo1Params,
              cryoMaterial,
              torsec->getRotation(),
              torsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );

    verbosityLevel &&
      G4cout << __func__ << " " << tssName << " Mass in kg: "
             << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg
             << G4endl;

    torsec = ts->getTSCryo<TorusSection>(TransportSolenoid::TSRegion::TS4,TransportSolenoid::TSRadialPart::OUT );
    std::array<double,5> ts4Cryo2Params { { torsec->rIn(), torsec->rOut(), torsec->torusRadius(), torsec->phiStart(), torsec->deltaPhi() } };
    tssName = "TS4OuterCryoShell";
    double ts4InsVacROut = torsec->rIn();
    nestTorus(tssName,
              ts4Cryo2Params,
              cryoMaterial,
              torsec->getRotation(),
              torsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );

    verbosityLevel &&
      G4cout << __func__ << " " << tssName << " Mass in kg: "
             << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg
             << G4endl;

    // Put in the insulating vacuum, which will serve as the mother volume
    // for the coils and coil assemblies (CAs).
    std::array<double,5> ts4CIVParams { { ts4InsVacRIn, ts4InsVacROut, torsec->torusRadius(), torsec->phiStart(), torsec->deltaPhi() } };

    tssName = "TS4CryoInsVac";
    nestTorus(tssName,
              ts4CIVParams,
              upstreamVacuumMaterial,
              torsec->getRotation(),
              torsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );


    if(tsThermalShield) {
      VolumeInfo useAsParent = _helper->locateVolInfo( "TS4CryoInsVac" );
      addThermalShield(*ts, useAsParent, config, TransportSolenoid::TSRegion::TS4,
                       thermalShieldMLIMaterial, thermalShieldMidMaterial, -1.);
    }


    // Build TS5.
    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS5,TransportSolenoid::TSRadialPart::IN );
    CLHEP::Hep3Vector globalVac5Position = CLHEP::Hep3Vector(strsec->getGlobal().x(),
                                                             strsec->getGlobal().y(),
                                                             strsec->getGlobal().z() );

    nestTubs( "TS5Vacuum",
              TubsParams( 0., ts->ts5InnerRadius(), strsec->getHalfLength() ),
              downstreamVacuumMaterial,
              strsec->getRotation(),
              globalVac5Position-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Green(),
              "TSCryo"
              );

    double ts5InsVacRIn = strsec->rOut();
    //account for possible overlap with coll53 (flange)
    double coll53HL = ts->getColl53().halfLength();
    double coll53RIn = ts->getColl53().rIn();
    double ts5InnerCryoHL = strsec->getHalfLength();
    double ts5InnerCryoOffset = 0.;
    if(coll53RIn < strsec->rOut()) {
      ts5InnerCryoHL -= coll53HL;
      ts5InnerCryoOffset = coll53HL;
    }
    tssName = "TS5InnerCryoShell";
    nestTubs( tssName,
              TubsParams( strsec->rIn(),
                          strsec->rOut(),
                          ts5InnerCryoHL ),
              cryoMaterial,
              strsec->getRotation(),
              strsec->getGlobal()-_hallOriginInMu2e - CLHEP::Hep3Vector(0.,0.,ts5InnerCryoOffset),
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );

    verbosityLevel &&
      G4cout << __func__ << " " << tssName << " Mass in kg: "
             << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg
             << G4endl;

    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS5,TransportSolenoid::TSRadialPart::OUT );
    tssName = "TS5OuterCryoShell";
    double ts5InsVacROut = strsec->rIn();
    nestTubs( tssName,
              TubsParams( strsec->rIn(),
                          strsec->rOut(),
                          strsec->getHalfLength() ),
              cryoMaterial,
              strsec->getRotation(),
              strsec->getGlobal()-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );

    verbosityLevel &&
      G4cout << __func__ << " " << tssName << " Mass in kg: "
             << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg
             << G4endl;

    if ( verbosityLevel > 0) {
      G4cout << __func__ << " TS5  OffsetInMu2e : " << strsec->getGlobal()   << G4endl;
      G4cout << __func__ << " TS5  rotation     : " << *(strsec->getRotation()) << G4endl;
    }


    // Put in the insulating vacuum, which will serve as the mother volume
    // for the coils and coil assemblies (CAs).

    double ts5InsVacHalfLen = strsec->getHalfLength() - ts->endWallD_halfLength();// - ts->getColl53().halfLength();
    double ts5InsVacOffset = ts->endWallD_halfLength();// + ts->getColl53().halfLength();
    CLHEP::Hep3Vector ts5CIVLoc = strsec->getGlobal()-_hallOriginInMu2e - CLHEP::Hep3Vector(0,0,ts5InsVacOffset);
    tssName =  "TS5CryoInsVac";
    nestTubs( tssName,
              TubsParams( ts5InsVacRIn,
                          ts5InsVacROut,
                          ts5InsVacHalfLen ),
              upstreamVacuumMaterial,
              strsec->getRotation(),
              ts5CIVLoc,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );


    if(tsThermalShield) {
      VolumeInfo useAsParent = _helper->locateVolInfo( "TS5CryoInsVac" );
      addThermalShield(*ts, useAsParent, config, TransportSolenoid::TSRegion::TS5,
                       thermalShieldMLIMaterial, thermalShieldMidMaterial, -1.);
    }

    // Build Rings (added April 1, 2015, David Norvil Brown)
    double rirs = ts->rInRingSide();
    double rors = ts->rOutRingSide();
    double trs = ts->thickRingSide();
    double rir = ts->rInRing();
    double ror = ts->rOutRing();
    double lr = ts->lengthRing();
    G4Material* ringMaterial = findMaterialOrThrow(ts->RingMaterial());
    std::vector<double> xr = ts->xRing();
    std::vector<double> yr = ts->yRing();
    std::vector<double> zr = ts->zRing();
    std::vector<double> thetasRing = ts->thetaRing();

    for ( unsigned int iRing = 0; iRing < xr.size(); iRing++ ) {
      // Let's make a mother volume first for each ring.
      std::ostringstream ringMotherName;
      ringMotherName << "TSRingMother" << iRing;
      CLHEP::HepRotation* ringRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
      double ringRotTheta = thetasRing[iRing]*CLHEP::degree;
      ringRotat->rotateY(ringRotTheta);
      CLHEP::HepRotation* noRotat (nullptr);

      double motherx = xr[iRing];
      double mothery = yr[iRing];
      double motherz = zr[iRing];

      VolumeInfo motherVol = nestTubs( ringMotherName.str(),
                                       TubsParams( rirs, rors, trs + lr/2.0 ),
                                       findMaterialOrThrow("G4_AIR"),
                                       ringRotat,
                                       CLHEP::Hep3Vector(motherx,mothery,motherz
                                                         ) - _hallOriginInMu2e,
                                       parent, 0, G4Color::Blue(),
                                       "TSCryo" );
      std::ostringstream leftName;
      leftName << "TSleftSideRing" << iRing;

      double lx = 0.0;
      double ly = 0.0;
      double lz = - lr/2.0 - trs/2.0;
      nestTubs( leftName.str(),
                TubsParams( rirs, rors, trs/2.0 ),
                ringMaterial,
                noRotat,
                CLHEP::Hep3Vector(lx,ly,lz),
                motherVol,
                0,
                G4Color::Blue(),
                "TSCryo"
                );

      std::ostringstream centerName;
      centerName << "TScenterRing" << iRing;

      nestTubs( centerName.str(),
                TubsParams( rir, ror, lr/2.0 ),
                ringMaterial,
                noRotat,
                CLHEP::Hep3Vector(0,0,0),
                motherVol,
                0,
                G4Color::Blue(),
                "TSCryo"
                );


      std::ostringstream rightName;
      rightName << "TSrightSideRing" << iRing;

      double rx = 0.0;
      double ry = 0.0;
      double rz = lr/2.0 + trs/2.0;
      if (iRing==0){

        G4VSolid *tsrightsolid= new G4Tubs( "TSrightSideRing0Full",rirs, rors, trs/2.0 ,0.,2*M_PI);
        AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
        // Get the beam pipe
        G4Tubs* beamPassTub = nullptr;
        CLHEP::HepRotation* turn = nullptr;
        CLHEP::Hep3Vector place;
        getBeamPipe(config, reg, beamPassTub, turn, place);

        CLHEP::Hep3Vector pos(rx,ry,rz);
        G4SubtractionSolid* arightSolid =
          new G4SubtractionSolid ( "TSrightSideRing0",
                                   tsrightsolid,
                                   beamPassTub,
                                   turn,
                                   place-pos-motherVol.centerInWorld);
        VolumeInfo tss(rightName.str(),
                       pos,
                       motherVol.centerInWorld);
        tss.solid = arightSolid;
        tss.solid->SetName(tss.name);


        finishNesting(tss,
                      ringMaterial,
                      0,
                      tss.centerInParent,
                      motherVol.logical,
                      0,
                      true,
                      G4Color::Blue(),
                      true,
                      true,
                      true,
                      false
                      );
      }
      else{
        nestTubs( rightName.str(),
                  TubsParams( rirs, rors, trs/2.0 ),
                  ringMaterial,
                  noRotat,
                  CLHEP::Hep3Vector(rx,ry,rz),
                  motherVol,
                  0,
                  G4Color::Blue(),
                  "TSCryo"
                  );
      }

    } // End of building rings

    // Build outer portion of support for TS3 antiproton stopping window
    // Will model for now as solid tubes

    int pbarAbsTS3Version = config.getInt("pbar.version",1);
    if ( pbarAbsTS3Version > 1 ) {

      std::ostringstream PabsSupOutName;
      PabsSupOutName << "PabsTS3SupOut";
      CLHEP::HepRotation* pasubRotat = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
      pasubRotat->rotateY(90.0*CLHEP::degree);
      double rPabsTS3SupOut_in  = rirs;
      double rPabsTS3SupOut_out = config.getDouble("pbar.support.outROut", rors);
      double rPabsTS3SupOut_HL  = config.getDouble("pbar.support.outHalfLength", trs);
      nestTubs( PabsSupOutName.str(),
                TubsParams( rPabsTS3SupOut_in, rPabsTS3SupOut_out, rPabsTS3SupOut_HL ),
                ringMaterial,
                pasubRotat,
                CLHEP::Hep3Vector(0,0,0)-_hallOriginInMu2e,
                parent,
                0,
                G4Color::Blue(),
                "TSCryo"
                );

      // Now do the next level in for the TS3 pabs window support - in cryo
      // (acts as endwall for TSu and TSd cryo sections.
      double rinner = config.getDouble("pbar.support.midRin")*CLHEP::mm;
      double router = config.getDouble("pbar.support.midRout")*CLHEP::mm;
      double halflen = config.getDouble("pbar.support.midThickness")*CLHEP::mm/2.0;

      Mu2eG4Helper* _helper = &(*art::ServiceHandle<Mu2eG4Helper>() );
      VolumeInfo useAsParent = _helper->locateVolInfo( "TS3CryoInsVac" );

      std::ostringstream PabsSupMidName;
      PabsSupMidName << "PabsTS3MidOut";
      nestTubs( PabsSupMidName.str(),
                TubsParams( rinner, router, halflen ),
                ringMaterial,
                0,
                CLHEP::Hep3Vector(0,0,0)- useAsParent.centerInMu2e(),
                useAsParent,
                0,
                G4Color::Blue(),
                "TSCryo"
                );

    } //end of " if pbarAbsTS3Version..."




    // Build downstream end wall of TS5
    // place flush with end of TS5
    CLHEP::Hep3Vector pos3( strsec->getGlobal().x(),
                            strsec->getGlobal().y(),
                            strsec->getGlobal().z()+strsec->getHalfLength()-ts->endWallD_halfLength() );

    tssName = "TS5DownstreamEndwall";
    nestTubs( tssName,
              TubsParams( ts->endWallD_rIn(),
                          ts->endWallD_rOut(),
                          ts->endWallD_halfLength() ),
              cryoMaterial,
              0,
              pos3-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
              "TSCryo"
              );
    //if defined a second component of the end wall
    if(ts->build_endWallD2()) {
      //place upstream of first component
      CLHEP::Hep3Vector pos32(pos3.x(), pos3.y(), pos3.z()-ts->endWallD_halfLength()-2.*ts->endWallD2_halfLength());
      VolumeInfo useAsParent = _helper->locateVolInfo( "TS5CryoInsVac" );
      pos32 -= useAsParent.centerInMu2e();
      tssName = "TS5DownstreamEndwall2";
      nestPolycone( tssName,
                    PolyconsParams(ts->endWallD2_z(),
                                   ts->endWallD2_rIn(),
                                   ts->endWallD2_rOut()),
                    cryoMaterial,
                    0,
                    pos32,
                    useAsParent,
                    0,
                    G4Color::Red(),
                    "TSCryo"
                    );
    }

    verbosityLevel &&
      G4cout << __func__ << " " << tssName << " Mass in kg: "
             << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg
             << G4endl;

    if ( verbosityLevel ) {
      G4cout << __func__ << " Downstream TS5 endwall at: " << pos3 << G4endl;
    }
  }

  //__________________________________
  //
  // CONSTRUCT Coil Assemblies (CA) (Modules), a torus approximation for now
  //__________________________________

  void constructCAs( VolumeInfo const& parent,
                     SimpleConfig const& config,
                     Beamline const& bl ) {

    typedef TransportSolenoid::TSCARegion::enum_type   tsCAReg_enum;

    const int  verbosityLevel = config.getInt("ts.cas.verbosityLevel", 0);

    TransportSolenoid const & ts       (bl.getTS());
    TorusSection      const * caTorsec (nullptr);
    ConeSection       const * caConsec (nullptr);
    StraightSection   const * caStrsec (nullptr);

    G4ThreeVector parentCenterInMu2e = parent.centerInMu2e();

    for ( unsigned iTS = tsCAReg_enum::TS1 ; iTS <= tsCAReg_enum::TS5 ; ++iTS ) {

      auto its = static_cast<TransportSolenoid::TSCARegion>(iTS);
      std::string const & caName  = its.name()+"CA";
      verbosityLevel && G4cout << __func__ << " constructing " << caName << G4endl;

      VolumeInfo useAsParent;
      Mu2eG4Helper* _helper = &(*art::ServiceHandle<Mu2eG4Helper>() );

      if ( its==tsCAReg_enum::TS2 || its==tsCAReg_enum::TS4 ) {

        // these sections are toruses
        // fixme, make the base class TSSection more general to be used here instead the two
        caTorsec = ts.getTSCA<TorusSection>(its);

        if ( its == tsCAReg_enum::TS2 ) {
          useAsParent = _helper->locateVolInfo( "TS2CryoInsVac" );
        } else if ( its == tsCAReg_enum::TS4 ) {
          useAsParent = _helper->locateVolInfo( "TS4CryoInsVac" );
        } else {  // Should never get here
          useAsParent = parent;
        }


        nestTorus(caName,
                  caTorsec->getParameters(),
                  findMaterialOrThrow(caTorsec->getMaterial()),
                  0,
                  caTorsec->getGlobal()-useAsParent.centerInMu2e(),
                  useAsParent,
                  0,
                  G4Color::Yellow(),
                  "TSCA"
                  );
        //                  caTorsec->getRotation(),
      } else if ( its==tsCAReg_enum::TS1 ){
        useAsParent = _helper->locateVolInfo( "TS1CryoInsVac" );
        caConsec = ts.getTSCA<ConeSection>(its);

        CLHEP::Hep3Vector V = caConsec->getGlobal() - useAsParent.centerInMu2e();
        CLHEP::HepRotation* rot = useAsParent.physical->GetRotation();
        V = (*rot) * V;


        G4VSolid *tssolid= new G4Cons("TS1CAFull",caConsec->rIn1(),caConsec->rOut1() ,caConsec->rIn2(),caConsec->rOut2(),caConsec->getHalfLength(),caConsec->phiStart(),caConsec->deltaPhi());

        AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
        // Get the beam pipe
        G4Tubs* beamPassTub = nullptr;
        CLHEP::HepRotation* turn = nullptr;
        CLHEP::Hep3Vector place;
        getBeamPipe(config, reg, beamPassTub, turn, place);

        CLHEP::Hep3Vector pos(caConsec->getGlobal());
        G4SubtractionSolid* aSolid =
          new G4SubtractionSolid (caName,
                                  tssolid,
                                  beamPassTub,
                                  turn,
                                  place-pos);

        G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();
        VolumeInfo tss(caName,
                       pos-_hallOriginInMu2e,
                       parent.centerInWorld);
        tss.solid = aSolid;
        tss.solid->SetName(tss.name);


        finishNesting(tss,
                      findMaterialOrThrow(caConsec->getMaterial()),
                      0,
                      V,
                      useAsParent.logical,
                      0,
                      true,
                      G4Color::Blue(),
                      true,
                      true,
                      true,
                      false
                      );

      } else if ( its==tsCAReg_enum::TS3u
                  || its==tsCAReg_enum::TS3d
                  || its==tsCAReg_enum::TS5
                  ) {

        if ( its == tsCAReg_enum::TS3u || its == tsCAReg_enum::TS3d ) {
          useAsParent = _helper->locateVolInfo( "TS3CryoInsVac" );
        } else if ( its == tsCAReg_enum::TS5 ) {
          useAsParent = _helper->locateVolInfo( "TS5CryoInsVac" );
        } else { // Should never get here
          useAsParent = parent;
        }

        caConsec = ts.getTSCA<ConeSection>(its);

        CLHEP::Hep3Vector V = caConsec->getGlobal() - useAsParent.centerInMu2e();
        CLHEP::HepRotation* rot = useAsParent.physical->GetRotation();
        V = (*rot) * V;


        nestCons(caName,
                 caConsec->getParameters(),
                 findMaterialOrThrow(caConsec->getMaterial()),
                 0,
                 V,
                 useAsParent,
                 0,
                 G4Color::Yellow(),
                 "TSCA"
                 );
        //                 caConsec->getRotation(),
      } else {

        useAsParent = _helper->locateVolInfo("TS3CryoInsVac");
        caStrsec = ts.getTSCA<StraightSection>(its);
        CLHEP::Hep3Vector V = caStrsec->getGlobal() - useAsParent.centerInMu2e();
        CLHEP::HepRotation* rot = useAsParent.physical->GetRotation();
        V = (*rot) * V;

        nestTubs(caName,
                 TubsParams(caStrsec->rIn(),
                            caStrsec->rOut(),
                            caStrsec->getHalfLength()),
                 findMaterialOrThrow(caStrsec->getMaterial()),
                 0,
                 V,
                 useAsParent,
                 0,
                 G4Color::Cyan(),
                 "TSCA"
                 );

        //caStrsec->getRotation(),

      }

      verbosityLevel
        && G4cout << __func__ << " " << caName << " Mass in kg: "
                  << _helper->locateVolInfo(caName).logical->GetMass()/CLHEP::kg
                  << G4endl;
    }

  }

  //__________________________________
  //
  // CONSTRUCT COILS
  //__________________________________

  void constructCoils( VolumeInfo const& parent,
                       SimpleConfig const& config,
                       Beamline const& bl ) {

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    const int  verbosityLevel      = config.getInt("ts.coils.verbosityLevel", 0);

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "TSCoils", "ts.coils"      );

    G4Material* coilMaterial = findMaterialOrThrow( bl.getTS().coil_material() );

    // Construct TS coils
    for ( unsigned iTS = TransportSolenoid::TSRegion::TS1 ; iTS <= TransportSolenoid::TSRegion::TS5 ; iTS++ ) {
      auto its = (TransportSolenoid::TSRegion::enum_type)iTS;
      int iC(0);
      for ( Coil const & coil : bl.getTS().getTSCoils( its ) ) {

        ostringstream coilname ; coilname << "TS" << iTS << "_Coil" << ++iC ;

        // Use insulating vacuums as mother volumes
        ostringstream ivName;
        ivName << "TS" << iTS << "CryoInsVac";

        Mu2eG4Helper* _helper = &(*art::ServiceHandle<Mu2eG4Helper>() );
        VolumeInfo useAsParent = _helper->locateVolInfo( ivName.str() );

        CLHEP::Hep3Vector V = coil.getGlobal() - useAsParent.centerInMu2e();


        CLHEP::HepRotation* rot = useAsParent.physical->GetRotation();
        V = (*rot) * V;
        CLHEP::HepRotation invrot = rot->inverse();

        if ( iTS == TransportSolenoid::TSRegion::TS2 ||
             iTS == TransportSolenoid::TSRegion::TS4 ) {

          CLHEP::HepRotation* twist = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
          twist->rotateX(90.0*CLHEP::degree);
          twist->rotateY(-coil.getRotation()->theta());


          nestTubs( coilname.str(),
                    TubsParams( coil.rIn(), coil.rOut(), coil.halfLength() ),
                    coilMaterial,
                    twist,
                    V,
                    useAsParent,
                    0,
                    G4Color::Green(),
                    "TSCoils"
                    );
        } else {

          nestTubs( coilname.str(),
                    TubsParams( coil.rIn(), coil.rOut(), coil.halfLength() ),
                    coilMaterial,
                    0,
                    V,
                    useAsParent,
                    0,
                    G4Color::Green(),
                    "TSCoils"
                    );
        }
        if ( verbosityLevel > 0 ) {
          G4cout << __func__ << " " << coilname.str() << " placed at: " << coil.getGlobal() << G4endl;
          G4cout << __func__ << "            rotation: " << -coil.getRotation()->getTheta()/CLHEP::degree << G4endl;
          G4cout << __func__ << "              params: " << coil.rIn() << " , " << coil.rOut() << " , " << 2*coil.halfLength() << G4endl;
        }

      }
    }

  }

  //__________________________________
  //
  // CONSTRUCT COLLIMATORS
  //__________________________________

  void constructCollimators( VolumeInfo const& parent,
                             SimpleConfig const& config,
                             Beamline const& bl ) {

    // Flags
    const int  verbosityLevel      = config.getInt("ts.coll.verbosityLevel", 0);

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "TSColl", "ts.coll"      );

    // Get collimators
    TransportSolenoid const& ts = bl.getTS();

    TSSection       const * ts1in = ts.getTSCryo(TransportSolenoid::TSRegion::TS1, TransportSolenoid::TSRadialPart::IN);
    TSSection       const * ts3in = ts.getTSCryo(TransportSolenoid::TSRegion::TS3, TransportSolenoid::TSRadialPart::IN);
    TorusSection    const * ts4in = ts.getTSCryo<TorusSection>(TransportSolenoid::TSRegion::TS4, TransportSolenoid::TSRadialPart::IN);
    StraightSection const * ts5in = ts.getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS5, TransportSolenoid::TSRadialPart::IN);

    CollimatorTS1 const& coll1   = ts.getColl1() ;
    CollimatorTS3 const& coll31  = ts.getColl31();
    CollimatorTS3 const& coll32  = ts.getColl32();
    CollimatorTS5 const& coll51  = ts.getColl51() ;
    CollimatorTS5 const& coll52  = ts.getColl52();
    CollimatorTS5 const& coll53  = ts.getColl53();

    // Get VDs
    GeomHandle<VirtualDetector> vdg;
    double vdHalfLength = vdg->getHalfLength()*CLHEP::mm;

    // Place collimator 1 (concentric cone which can be a cylinder when r1==r2
    // and a cylinder placed in TS1Vacuum)

    // the cone (which can be a tube/cylinder) inside the outer tube/cylinder
    double coll11rout = coll1.rIn3();
    if(coll1.rOut1() > 0.) //allow the inner liner to have a different outer radius than the next layer's inner radius
      coll11rout = coll1.rOut1();

    double coll1Param1[7] = { coll1.rIn1(), coll11rout,
      coll1.rIn2(), coll11rout,
      coll1.halfLength() - 2.*vdHalfLength,
      0.0, CLHEP::twopi };

    Mu2eG4Helper* _helper = &(*art::ServiceHandle<Mu2eG4Helper>() );

    CLHEP::Hep3Vector parentPosW    = _helper->locateVolInfo("TS1Vacuum").centerInWorld;
    CLHEP::Hep3Vector parentPosM    = _helper->locateVolInfo("TS1Vacuum").centerInMu2e();


    if ( verbosityLevel > 0 ) {
      G4cout << __func__ << " Coll1 OffsetInW     : " << coll1.getLocal() + parentPosW << G4endl;
      G4cout << __func__ << " Coll1 OffsetInMu2e  : " << coll1.getLocal() + parentPosM << G4endl;
      G4cout << __func__ << " Coll1 Extent        :[ " << coll1.getLocal().z() - coll1.halfLength() + parentPosM.z() <<","
             << coll1.getLocal().z() + coll1.halfLength()  + parentPosM.z() << "]" << G4endl;
    }



    nestCons( "Coll11",
              coll1Param1,
              findMaterialOrThrow( coll1.material1() ),
              0,
              coll1.getLocal(),
              _helper->locateVolInfo("TS1Vacuum"),
              0,
              G4Color::Cyan(),
              "TSColl"
              );

    double tmpRout = coll1.rOut();
    if ( coll1.rIn4() > 1.0 ) tmpRout = coll1.rIn4();

    TubsParams coll1Param2 ( coll1.rIn3(),  tmpRout, coll1.halfLength()-2.*vdHalfLength);

    nestTubs( "Coll12",
              coll1Param2,
              findMaterialOrThrow( coll1.material2() ),
              0,
              coll1.getLocal(),
              _helper->locateVolInfo("TS1Vacuum"),
              0,
              G4Color::Blue(),
              "TSColl"
              );

    if ( coll1.rIn4() > 1.0 && coll1.rOut() > coll1.rIn4() ) {
      // Make the sheath
      TubsParams coll1Param3 ( coll1.rIn4(),  coll1.rOut(), coll1.halfLength()-2.*vdHalfLength);
      nestTubs( "Coll13",
                coll1Param3,
                findMaterialOrThrow( coll1.material3() ),
                0,
                coll1.getLocal(),
                _helper->locateVolInfo("TS1Vacuum"),
                0,
                G4Color::Blue(),
                "TSColl"
                );

    } // end of adding sheath to Coll1


    if ( verbosityLevel > 0) {
      G4cout << __func__ << " TS1  OffsetInMu2e    : " << ts1in->getGlobal()       << G4endl;
      G4cout << __func__ << " Coll1 local offset   : " << ts.getColl1().getLocal() << G4endl;
      G4cout << __func__ << " TS1  Rotation        : " << ts1in->getRotation()     << G4endl;
    }

    // Place collimator 3

    // Collimator 3 has peculiar shape, described in doc_db 853.
    // Construct this shape using boolean functions on solids

    // First, construct hole; make it slightly longer than any collimator
    double hDz = coll31.halfLength();
    if( hDz<coll32.halfLength() ) hDz=coll32.halfLength();
    // Hole is the intersection of box and tube
    G4Box* coll3_hole_box = new G4Box("coll3_hole_box",
                                      coll31.holeRadius()+5.0,coll31.holeHalfHeight(),hDz+1.0);
    // make the tube longer than the box to avoid overlapping surfaces
    G4Tubs* coll3_hole_circle = new G4Tubs("coll3_hole_circle",
                                           0.0,coll31.holeRadius(),hDz+2.0,
                                           0.0, CLHEP::twopi );
    G4IntersectionSolid* coll3_hole = new G4IntersectionSolid("coll3_hole",
                                                              coll3_hole_box,
                                                              coll3_hole_circle);


    // Make collimators themselves. At this moment the collimators
    // coll31 and coll32 are the same size. But it is possible to make them
    // different length. Therefore two solids are created, but the same hole
    // is subtracted from both solids.

    VolumeInfo coll31Info;
    VolumeInfo coll32Info;

    coll31Info.name = "Coll31";
    coll32Info.name = "Coll32";

    if ( verbosityLevel > 0 ) {
      CLHEP::Hep3Vector parentPosW    = _helper->locateVolInfo("TS3Vacuum").centerInWorld;
      CLHEP::Hep3Vector parentPosM    = _helper->locateVolInfo("TS3Vacuum").centerInMu2e();
      G4cout << __func__ << " Coll31 OffsetInW     : " << coll31.getLocal() + parentPosW << G4endl;
      G4cout << __func__ << " Coll31 OffsetInMu2e  : " << coll31.getLocal() + parentPosM << G4endl;
      G4cout << __func__ << " Coll31 Extent        :[ " << coll31.getLocal().z() - coll31.halfLength() + parentPosM.z() <<","
             << coll31.getLocal().z() + coll31.halfLength()  + parentPosM.z() << "]" << G4endl;
      G4cout << __func__ << " ts innerRadius      : " << ts.ts3InnerRadius() << G4endl;
      G4cout << __func__ << " coll31 outerRadius  : " << coll31.rOut()    << G4endl;
    }

    G4Tubs* coll31_mother = new G4Tubs("Coll31_mother",
                                       0, coll31.rOut(), coll31.halfLength()-2.*vdHalfLength,
                                       0.0, CLHEP::twopi );

    G4Tubs* coll32_mother = new G4Tubs("Coll32_mother",
                                       0, coll32.rOut(), coll32.halfLength()-2.*vdHalfLength,
                                       0.0, CLHEP::twopi );

    coll31Info.solid = new G4SubtractionSolid(coll31Info.name,
                                              coll31_mother,
                                              coll3_hole,
                                              0,
                                              G4ThreeVector(0,coll31.holeDisplacement(),0));

    coll32Info.solid = new G4SubtractionSolid(coll32Info.name,
                                              coll32_mother,
                                              coll3_hole,
                                              0,
                                              G4ThreeVector(0,coll32.holeDisplacement(),0));

    // Now use finishNesting to place collimators 31 and 32
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    G4RotationMatrix* coll31Rot = reg.add(G4RotationMatrix());
    G4RotationMatrix* coll32Rot = reg.add(G4RotationMatrix());
    coll31Rot->rotateZ(coll31.rotationAngle()*CLHEP::degree);
    coll32Rot->rotateZ(coll32.rotationAngle()*CLHEP::degree);

    finishNesting(coll31Info,
                  findMaterialOrThrow( coll31.material() ),
                  coll31Rot,
                  coll31.getLocal(),
                  _helper->locateVolInfo("TS3Vacuum").logical,
                  0,
                  G4Color::Gray(),
                  "TSColl");

    finishNesting(coll32Info,
                  findMaterialOrThrow( coll32.material() ),
                  coll32Rot,
                  coll32.getLocal(),
                  _helper->locateVolInfo("TS3Vacuum").logical,
                  0,
                  G4Color::Gray(),
                  "TSColl");

    // **************************************************************
    // Now place the "flashblocks" for tests of mitigating beam flash
    // **************************************************************

    if (coll31.useFlashBlock()) {
      std::vector<double> boxPars = { coll31.flashBlockWidth()/2.0*CLHEP::mm,
        coll31.flashBlockHeight()/2.0*CLHEP::mm,
        coll31.flashBlockLength()/2.0*CLHEP::mm};

      CLHEP::Hep3Vector displaceFB(coll31.flashBlockTranOff()*CLHEP::mm,
                                   coll31.holeDisplacement() - coll31.holeHalfHeight() + coll31.flashBlockHeight()/2.0*CLHEP::mm,
                                   coll31.flashBlockLength()/2.0 - coll31.halfLength());

      nestBox( "flashBlockUp",
               boxPars,
               findMaterialOrThrow(coll31.flashBlockMaterial()),
               coll31Rot,
               coll31.getLocal()+displaceFB,
               _helper->locateVolInfo("TS3Vacuum").logical,
               0,
               G4Colour::Gray(),
               "TSColl");
    }

    if (coll32.useFlashBlock()) {
      std::vector<double> boxPars = { coll32.flashBlockWidth()/2.0*CLHEP::mm,
        coll32.flashBlockHeight()/2.0*CLHEP::mm,
        coll32.flashBlockLength()/2.0*CLHEP::mm};

      CLHEP::Hep3Vector displaceFB(coll32.flashBlockTranOff()*CLHEP::mm,
                                   coll32.holeDisplacement() - coll32.holeHalfHeight() + coll32.flashBlockHeight()/2.0*CLHEP::mm,
                                   coll32.flashBlockLength()/2*CLHEP::mm - coll32.halfLength());

      nestBox( "flashBlockDn",
               boxPars,
               findMaterialOrThrow(coll32.flashBlockMaterial()),
               coll32Rot,
               coll32.getLocal()+displaceFB,
               _helper->locateVolInfo("TS3Vacuum").logical,
               0,
               G4Colour::Gray(),
               "TSColl");
    }


    // Now add a Recorder at the Coll31 exit and Coll32 entrance
    // (do not use VirtualDetector because of _ in its volume name)
    TubsParams coll31OutRecordParam ( 0,  ts.ts3InnerRadius(), vdHalfLength );
    G4ThreeVector coll31OutRecordTrans( coll31.getLocal().x(),
                                        coll31.getLocal().y(),
                                        coll31.getLocal().z() + coll31.halfLength() + vdHalfLength );
    nestTubs( "Coll31OutRecord",
              coll31OutRecordParam,
              findMaterialOrThrow(ts.downstreamVacuumMaterial()),
              coll31Rot,
              coll31OutRecordTrans,
              _helper->locateVolInfo("TS3Vacuum"),
              0,
              G4Color::Blue(),
              "TSColl"
              );

    TubsParams coll32InRecordParam ( 0,  coll32.rOut(), vdHalfLength );
    G4ThreeVector coll32InRecordTrans( coll32.getLocal().x(),
                                       coll32.getLocal().y(),
                                       coll32.getLocal().z() - coll32.halfLength() - vdHalfLength );
    nestTubs( "Coll32InRecord",
              coll32InRecordParam,
              findMaterialOrThrow(ts.downstreamVacuumMaterial()),
              coll32Rot,
              coll32InRecordTrans,
              _helper->locateVolInfo("TS3Vacuum"),
              0,
              G4Color::Blue(),
              "TSColl"
              );

    if ( verbosityLevel > 0) {
      G4cout << __func__ << " TS3  OffsetInMu2e   : " << ts3in->getGlobal() << G4endl;
      G4cout << __func__ << " Coll31 local offest : " << coll31.getLocal() << G4endl;
      G4cout << __func__ << " Coll32 local offset : " << coll32.getLocal() << G4endl;
      G4cout << __func__ << " TS3  Rotation       : " << ts3in->getRotation() << G4endl;
    }

    // Place collimator 5

    if ( verbosityLevel > 0) {
      G4cout << __func__ << " TS5  OffsetInMu2e  : " << ts5in->getGlobal()   << G4endl;
      G4cout << __func__ << " Coll5 local offset : " << coll51.getLocal()    << G4endl;
      G4cout << __func__ << " TS5  Rotation      : " << ts5in->getRotation() << G4endl;
    }

    CLHEP::Hep3Vector coll5OffsetInMu2e = ts5in->getGlobal() +
      ( ( ts5in->getRotation() != 0x0 ) ?
        *(ts5in->getRotation()) * coll51.getLocal() :
        coll51.getLocal() );

    if ( verbosityLevel > 0) {
      G4cout << __func__ << "  coll5OffsetInMu2e    : "    << coll5OffsetInMu2e << G4endl;
      G4cout << __func__ << "  Coll5 calc local offset : " << coll5OffsetInMu2e - ts5in->getGlobal() << G4endl;
    }

    // the most outer part (with Virtual Detectors on the outer surfaces of the Coll5)
    TubsParams coll5Param1 ( coll51.rIn(),  coll51.rOut() - 2.*vdHalfLength, coll51.halfLength()-2.*vdHalfLength);
    nestTubs( "Coll51",
              coll5Param1,
              findMaterialOrThrow( coll51.material() ),
              0,
              coll51.getLocal(),
              _helper->locateVolInfo("TS5Vacuum"),
              0,
              G4Color::Blue(),
              "TSColl"
              );


    double coll52HL = coll52.halfLength()-2.*vdHalfLength;
    double coll52HLInTS4 = -1.;
    auto coll52Origin = coll52.getLocal();
    if(coll52Origin.z() - coll52HL < -1.*ts5in->getHalfLength()) { //extends into TS4
      coll52HLInTS4 = coll52HL;
      double lengthOutTS5 = coll52HL - ts5in->getHalfLength() - coll52Origin.z(); //amount in TS4
      coll52HL -= lengthOutTS5/2.; //fill to end of upstream end of TS5
      coll52Origin.setZ(coll52Origin.z() + lengthOutTS5/2.); //offset the TS5 section to compensate
      coll52HLInTS4 = lengthOutTS5/2.; //left over HL in TS4
    }

    TubsParams coll5Param2 ( coll52.rIn(),  coll52.rOut(), coll52HL);

    nestTubs( "Coll52",
              coll5Param2,
              findMaterialOrThrow( coll52.material() ),
              0,
              coll52Origin,
              _helper->locateVolInfo("TS5Vacuum"),
              0,
              G4Color::Blue(),
              "TSColl"
              );
    if(coll52HLInTS4 > 0.) {
      TubsParams coll5Param2_TS4 ( coll52.rIn(),  coll52.rOut(), coll52HLInTS4);
      CLHEP::HepRotation* rotColl52 = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
      rotColl52->rotateX(90.*CLHEP::degree);
      nestTubs( "Coll52InTS4",
                coll5Param2_TS4,
                findMaterialOrThrow( coll52.material() ),
                rotColl52,
                CLHEP::Hep3Vector(-ts4in->torusRadius(),coll52HLInTS4,0.), //downstream end of TS4
                _helper->locateVolInfo("TS4Vacuum"),
                0,
                G4Color::Blue(),
                "TSColl"
                );
    }

    int coll53version = coll53.version();
    TubsParams coll5Param3 ( coll53.rIn(),  coll53.rOut(), coll53.halfLength()-2.*vdHalfLength);
    VolumeInfo useAsParent;
    G4ThreeVector _useAsParentOriginInMu2e;
    CLHEP::Hep3Vector offset(0.,0.,0.);
    StraightSection   const * strsec (nullptr);
    strsec = ts.getTSVacuum<StraightSection>(TransportSolenoid::TSRegion::TS5);
    if(coll53version == 1) {
      // Build Coll53 flange outside the TS5Vacuum
      useAsParent = _helper->locateVolInfo("TS5CryoInsVac");
      _useAsParentOriginInMu2e = useAsParent.centerInMu2e();
    } else if(coll53version > 1) { //place in hall air, flush to TS end
      //adjust origin from inside TS5 to end of it
      offset.setZ(2.*(ts.endWallD_halfLength()));
      if(ts.ts5InnerRadius() > coll53.rIn()) { //split between ts5 vacuum and hall air
        useAsParent = _helper->locateVolInfo("TS5Vacuum");
        _useAsParentOriginInMu2e = useAsParent.centerInMu2e();
        coll5Param3  = TubsParams(coll53.rIn(), ts.ts5InnerRadius(), coll53.halfLength()-2.*vdHalfLength);
        nestTubs( "Coll53InTS5",
                  coll5Param3,
                  findMaterialOrThrow( coll53.material() ),
                  0,
                  strsec->getGlobal()-_useAsParentOriginInMu2e+coll53.getLocal() + offset,
                  useAsParent,
                  0,
                  G4Color::Blue(),
                  "TSColl"
                  );
        coll5Param3  = TubsParams(ts.ts5InnerRadius(), coll53.rOut(), coll53.halfLength()-2.*vdHalfLength);
      }
      useAsParent = parent; //the rest of the flange is in hall air
      _useAsParentOriginInMu2e = parent.centerInMu2e();
    }
    nestTubs( "Coll53",
              coll5Param3,
              findMaterialOrThrow( coll53.material() ),
              0,
              strsec->getGlobal()-_useAsParentOriginInMu2e+coll53.getLocal() + offset,
              useAsParent,
              0,
              G4Color::Blue(),
              "TSColl"
              );
  }

  //__________________________________
  //
  // CONSTRUCT DEGRADER
  //__________________________________

  void constructDegrader( VolumeInfo const& parent,
                          SimpleConfig const & config,
                          Beamline const& bl ) {

    // Add muon degrader to the center of TS5
    // Degrader is made of several layers.
    // Each layer is intersection of cylinder and trapezoid.

    // Flag
    int const verbosityLevel = config.getInt("muondegrader.verbosityLevel", 0);

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "MuonDegrader", "muondegrader" );

    bool addDegrader  = config.getBool("muondegrader.build",false);
    vector<double> degraderR, degraderDZB, degraderDZT, degraderPhi;
    config.getVectorDouble("muondegrader.R", degraderR, vector<double>() );
    config.getVectorDouble("muondegrader.DZB", degraderDZB, vector<double>() );
    config.getVectorDouble("muondegrader.DZT", degraderDZT, vector<double>() );
    config.getVectorDouble("muondegrader.Phi", degraderPhi, vector<double>() );

    CollimatorTS5 const& coll5 = bl.getTS().getColl51();

    if( degraderR.size()!=degraderDZB.size() || degraderR.size()!=degraderDZT.size() ||
        degraderR.size()!=degraderPhi.size() ) {
      G4cout << __func__ << " Warning: MuonDegrader is not build - dimensions don't match." << G4endl;
      addDegrader = false;
    }

    if( addDegrader ) {

      G4Material* degraderMaterial  = findMaterialOrThrow( config.getString("muondegrader.materialName") );

      for( unsigned int i=0; i<degraderR.size(); ++i ) {

        VolumeInfo degraderInfo;

        ostringstream dname;  dname << "MuonDegrader" << i;
        degraderInfo.name = dname.str();

        ostringstream dsname1;  dsname1 << "MuonDegrader_disk" << i;
        ostringstream dsname2;  dsname2 << "MuonDegrader_trd" << i;

        double R1 = degraderR[i];
        double R2 = ( i<=0 ) ? coll5.rIn() : degraderR[i-1];
        G4Tubs *degrader_disk = new G4Tubs(dsname1.str(),R1,R2,
                                           degraderDZB[i]/2.0,
                                           (1.5-degraderPhi[i])*CLHEP::pi,
                                           degraderPhi[i]*CLHEP::twopi);

        G4Trd *degrader_trd = new G4Trd(dsname2.str(),
                                        coll5.rIn(),coll5.rIn(),
                                        degraderDZT[i]/2.0,degraderDZB[i]/2.0,
                                        R2/2);

        AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
        G4RotationMatrix* degRot = reg.add(G4RotationMatrix());
        degRot->rotateX(-90.0*CLHEP::degree);
        G4ThreeVector degTrans(0.0,-R2/2,0.0);

        degraderInfo.solid = new G4IntersectionSolid(degraderInfo.name,
                                                     degrader_disk,
                                                     degrader_trd,
                                                     degRot,
                                                     degTrans);

        finishNesting(degraderInfo,
                      degraderMaterial,
                      0,
                      coll5.getLocal(),
                      parent.logical,
                      0,
                      G4Color::Blue()
                      );

        if ( verbosityLevel > 0 ) {
          G4cout << __func__ << " Degrader constructed at: " << coll5.getLocal() << " wrt. TS5 " << G4endl;
        }

      }

    }

  }

  //__________________________________
  //
  // CONSTRUCT PBAR WINDOW
  //__________________________________

  void constructPbarWindow( VolumeInfo const& parent,
                            SimpleConfig const& config,
                            Beamline const& bl ) {
    // ******* These are notes for version 1 of the pbar window ********
    // Place Pbar absorber between Coll31 and Coll32
    // Pbar absorber is made of two pieces:
    //  -- vacuum wall, which covers the whole inner part of TS3
    //     it is controlled by pbar.* parameters
    //  -- wedge, which starts near center and extends upward
    //     it is controlled by pbarwedge.* parameters
    // ******* In version two, there are changes *****
    // - support structure is ~1 cm thick stainless with a window shaped
    //   like that of the COL3u and Col3d windows.
    // - wedge is shaped like the
    //   hole in the support structure.
    // ******* In version three, there are some more changes
    // - The window in the support structure is circular with the circular
    //   disk to be placed inside, with a radius less than the inner cryo
    //   shell.
    // - The wedge has a rectangular projection in the y-z plane.
    //   It is made of a series of rectangular sheets forming a
    //   "stairstep" structure.
    // ******* In Version 4, do the same as version 3 except allow the steps to have different
    //   thicknesses.  This will allow us to implement the geometry of doc-db 17519 p. 34

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "PbarAbs", "pbar" );

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    // Throw exception if pbarwedge.build is used - way out of date!
    if ( config.hasName("pbarwedge.build") )
      {
        throw cet::exception("GEOM")<<
          " Variable pbarwedge.build is now deprecated. \n" <<
          " To use pbar wedge specify: pbar.Type = \"wedge\" \n" ;
      }

    PbarWindow const & pbarWindow = bl.getTS().getPbarWindow();
    G4Material* pbarMaterial  = findMaterialOrThrow( pbarWindow.material() );

    // First, ascertain which version this is, along with other config info
    int pbarAbsTS3Version = pbarWindow.version();
    int const verbosityLevel = config.getInt("pbar.verbosityLevel", 0);

    if (verbosityLevel > 0){
      G4cout << "pbarWindow.shape() = " << pbarWindow.shape() << G4endl;
    }

    if ( pbarAbsTS3Version == 1 ) {
      // -- vacuum wall

      G4cout << "inside version 1 " << G4endl;
      if (verbosityLevel > 0) G4cout << "TS3 pbar windows HalfLength : " << pbarWindow.halfLength() << G4endl;

      if ( pbarWindow.shape() == "wedge" ||
           pbarWindow.shape() == "disk" ) {

        double pbarParams[5]  = { 0.0, pbarWindow.rOut(), pbarWindow.halfLength(), 0.0, CLHEP::twopi };

        nestTubs( "PbarAbs",
                  pbarParams,
                  pbarMaterial,
                  0,
                  pbarWindow.getLocal(),
                  parent,
                  0,
                  G4Color::Yellow()
                  );

      }

      //      if( pbarWindow.shape() == "wedge" )
      if ( pbarWindow.shape() != "disk")
        {
          if( pbarWindow.shape() == "wedge")
            {
              // -- pbar wedge
              double pbarWedge_y0  = pbarWindow.getY0();
              double pbarWedge_y1  = pbarWindow.getY1();
              double pbarWedge_dz0 = pbarWindow.getDZ0();
              double pbarWedge_dz1 = pbarWindow.getDZ1();

              VolumeInfo pbarWedgeInfo;

              pbarWedgeInfo.name = "PbarAbsWedge";

              double pbarWedge_dz = ( pbarWedge_dz0<pbarWedge_dz1 ) ? pbarWedge_dz1 : pbarWedge_dz0;
              double pbarWedge_h = pbarWedge_y1 - pbarWedge_y0;

              double pbarWedge_dy = (pbarWedge_y1 + pbarWedge_y0)/2.;

              G4Tubs *pbarWedge_disk = new G4Tubs("PbarAbsWedge_disk",
                                                  0,bl.getTS().innerRadius(),pbarWedge_dz/2.,0,CLHEP::twopi);

              G4Trd *pbarWedge_trd = new G4Trd("PbarAbsWedge_trd",
                                               bl.getTS().innerRadius(),bl.getTS().innerRadius(),
                                               pbarWedge_dz0/2.,pbarWedge_dz1/2.,
                                               pbarWedge_h/2.);

              AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
              G4RotationMatrix* pbarWedgeRot = reg.add(G4RotationMatrix());
              pbarWedgeRot->rotateX(90.0*CLHEP::degree);
              G4ThreeVector pbarWedgeTrans(0.0,pbarWedge_dy,0.0);

              pbarWedgeInfo.solid = new G4IntersectionSolid(pbarWedgeInfo.name,
                                                            pbarWedge_disk,
                                                            pbarWedge_trd,
                                                            pbarWedgeRot,
                                                            pbarWedgeTrans);

              finishNesting(pbarWedgeInfo,
                            pbarMaterial,
                            0,
                            G4ThreeVector(0.,0.,pbarWedge_dz/2+pbarWindow.halfLength() ),
                            parent.logical,
                            0,
                            G4Color::Yellow(),
                            "PbarAbs"
                            );

            }
          else if ( pbarWindow.shape() == "polycone" )
            {
              // Define polycone parameters
              vector<double> tmp_zPlanesDs3 {-0.50,-0.06,0.06,0.50};
              vector<double> tmp_rOuterDs3  (4,239.5);
              vector<double> tmp_rInnerDs3  {239.5,0.,0.,239.5};


              CLHEP::Hep3Vector polyPositionInMu2e = bl.getTS().getTSCryo(TransportSolenoid::TSRegion::TS3,
                                                                          TransportSolenoid::TSRadialPart::IN)->getGlobal();

              nestPolycone( "PbarAbsPolycone",
                            PolyconsParams(tmp_zPlanesDs3,
                                           tmp_rInnerDs3,
                                           tmp_rOuterDs3 ),
                            pbarMaterial,
                            0,
                            polyPositionInMu2e - parent.centerInMu2e(),
                            parent,
                            0,
                            G4Colour::Yellow(),
                            "PbarAbs"
                            );
            }
          else if ( pbarWindow.shape() != "disk" )
            {
              throw cet::exception("GEOM")<<
                " Incorrect pbar window geometry requested! \n << pbarWindow.shape() = in version 1" << pbarWindow.shape() <<"\n";
            }
        }
    }  // end of if ( pbarAbsTS3Version == 1 )
    else if ( pbarAbsTS3Version == 2 ) {
      // =============== Now Version 2 of pbarAbs in TS3! ==============
      // First, build the subtraction shape to represent the hole in the
      // support.  Based on code in Collimator build function
      // Get collimators

      G4cout << "inside version 2" << G4endl;
      TransportSolenoid const& ts = bl.getTS();
      CollimatorTS3 const& coll31  = ts.getColl31();
      // First, construct hole; make it slightly longer than the support
      double hDz = config.getDouble("pbar.support.innerHalflength")* CLHEP::mm;

      VolumeInfo supportInfo;
      supportInfo.name = "pBarAbsSupport";

      // Hole is the intersection of box and tube
      G4Box* support_hole_box = new G4Box("support_hole_box",
                                          coll31.holeRadius()+5.0,coll31.holeHalfHeight(),hDz+1.0);
      // make the tube longer than the box to avoid overlapping surfaces
      G4Tubs* support_hole_circle = new G4Tubs("support_hole_circle",
                                               0.0,coll31.holeRadius(),hDz+2.0,
                                               0.0, CLHEP::twopi );
      G4IntersectionSolid* support_hole = new G4IntersectionSolid("support_hole",
                                                                  support_hole_box,
                                                                  support_hole_circle);

      // Now make the actual support
      G4Tubs* support_mother = new G4Tubs("PbarSupport_mother",
                                          0, coll31.rOut(), hDz,
                                          0.0, CLHEP::twopi );

      supportInfo.solid = new G4SubtractionSolid("pBarTS3Support",
                                                 support_mother,
                                                 support_hole,
                                                 0,
                                                 G4ThreeVector(0,
                                                               config.getDouble("pbar.support.holeDisp"),0));

      CLHEP::HepRotation * supportRot (nullptr);
      //      supportRot->rotateY(90.0*CLHEP::degree);

      finishNesting(supportInfo,
                    findMaterialOrThrow(config.getString("pbar.support.material")),
                    supportRot,
                    G4ThreeVector(0,0,0),
                    parent.logical,
                    0,
                    G4Color::Gray(),
                    "PbarAbs");

      // -- vacuum wall

      if (verbosityLevel > 0) G4cout << "TS3 pbar windows HalfLength : " << pbarWindow.halfLength() << G4endl;
      if ( pbarWindow.shape() == "wedge" ||
           pbarWindow.shape() == "disk" ) {

        VolumeInfo pbarDiskInfo;
        pbarDiskInfo.name = "PbarAbsDisk";

        // Helper info
        double pbarWedge_y0  = pbarWindow.getY0();
        double pbarWedge_y1  = pbarWindow.getY1();
        double pbarWedge_dy = (pbarWedge_y1 + pbarWedge_y0)/2.;

        G4Tubs *pbarAbs_disk = new G4Tubs("PbarAbs_disk",
                                          0.0 ,pbarWindow.rOut(),
                                          pbarWindow.halfLength(),
                                          0.0,CLHEP::twopi);

        pbarDiskInfo.solid = new G4IntersectionSolid(pbarDiskInfo.name,
                                                     support_hole,
                                                     pbarAbs_disk,
                                                     0,
                                                     G4ThreeVector(0,0,0));

        finishNesting(pbarDiskInfo,
                      pbarMaterial,
                      0,
                      G4ThreeVector(0,pbarWedge_dy,0),
                      parent.logical,
                      0,
                      G4Color::Yellow(),
                      "PbarAbs"
                      );
      }

      if( pbarWindow.shape() == "wedge" )
        {
          // -- pbar wedge
          double pbarWedge_y0  = pbarWindow.getY0();
          double pbarWedge_y1  = pbarWindow.getY1();
          double pbarWedge_dz0 = pbarWindow.getDZ0();
          double pbarWedge_dz1 = pbarWindow.getDZ1();
          double pbarWedge_dz = ( pbarWedge_dz0<pbarWedge_dz1 ) ? pbarWedge_dz1 : pbarWedge_dz0;
          double pbarWedge_offsetZ = pbarWindow.getWedgeZOffset() + pbarWedge_dz/2.0;

          VolumeInfo pbarWedgeInfo;

          pbarWedgeInfo.name = "PbarAbsWedge";

          double pbarWedge_h = pbarWedge_y1 - pbarWedge_y0;

          double pbarWedge_dy = (pbarWedge_y1 + pbarWedge_y0)/2.;

          G4Trd *pbarWedge_trd = new G4Trd("PbarAbsWedge_trd",
                                           bl.getTS().innerRadius(),bl.getTS().innerRadius(),
                                           pbarWedge_dz0/2.,pbarWedge_dz1/2.,
                                           pbarWedge_h/2.);

          AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
          G4RotationMatrix* pbarWedgeRot = reg.add(G4RotationMatrix());
          pbarWedgeRot->rotateX(90.0*CLHEP::degree);
          G4ThreeVector pbarWedgeTrans(0.0,pbarWedge_dy,pbarWedge_offsetZ);

          pbarWedgeInfo.solid = new G4IntersectionSolid(pbarWedgeInfo.name,
                                                        support_hole,
                                                        pbarWedge_trd,
                                                        pbarWedgeRot,
                                                        G4ThreeVector(0,0,0));

          finishNesting(pbarWedgeInfo,
                        pbarMaterial,
                        0,
                        pbarWedgeTrans,
                        parent.logical,
                        0,
                        G4Color::Yellow(),
                        "PbarAbs"
                        );
        } //end of if ( pbarWindow.shape == wedge )
    }  else if ( pbarAbsTS3Version == 3 || pbarAbsTS3Version == 4) {

      if ( verbosityLevel > 2) {
        G4cout << __func__ <<"inside version " << pbarAbsTS3Version << G4endl;
      }
      // =============== Now Version 3 of pbarAbs in TS3! ==============
      // Get collimators (we use the coll31 info)
      CollimatorTS3 const& coll31  = bl.getTS().getColl31();
      // First, construct hole; make it slightly longer than the support
      double hDz = config.getDouble("pbar.support.innerHalflength")* CLHEP::mm;

      VolumeInfo supportInfo;
      supportInfo.name = "pBarAbsSupport";


      // Now make the actual support
      G4Tubs* support_mother = new G4Tubs("PbarSupport_mother",
                                          pbarWindow.diskRadius(),
                                          coll31.rOut(), hDz,
                                          0.0, CLHEP::twopi );

      supportInfo.solid = support_mother;

      CLHEP::HepRotation * supportRot (nullptr);

      finishNesting(supportInfo,
                    findMaterialOrThrow(config.getString("pbar.support.material")),
                    supportRot,
                    G4ThreeVector(0,0,0),
                    parent.logical,
                    0,
                    G4Color::Gray(),
                    "PbarAbs");

      // -- vacuum wall


      if ( pbarWindow.shape() == "wedge" ||
           pbarWindow.shape() == "disk" ) {

        VolumeInfo pbarDiskInfo;
        pbarDiskInfo.name = "PbarAbsDisk";

        if (verbosityLevel > 0) {
          G4cout << "TS3 pbar window thickness : " << pbarWindow.halfLength()*2. << G4endl;
          if (verbosityLevel > 1){
            G4cout << " inside wedge or disk" << G4endl;
          }
        }

        pbarDiskInfo.solid = new G4Tubs("PbarAbs_disk",
                                        0.0 ,pbarWindow.diskRadius(),
                                        pbarWindow.halfLength(),
                                        0.0,CLHEP::twopi);

        finishNesting(pbarDiskInfo,
                      pbarMaterial,
                      0,
                      G4ThreeVector(0,0,0),
                      parent.logical,
                      0,
                      G4Color::Yellow(),
                      "PbarAbs"
                      );

      }

      if( pbarWindow.shape() == "wedge" )
        {
          // Helper info
          // -- pbar wedge
          double pbarWedge_y1  = pbarWindow.getY1();
          double pbarWedge_offsetZ = pbarWindow.getWedgeZOffset();
          VolumeInfo pbarWedgeInfo;

          pbarWedgeInfo.name = "PbarAbsWedge";

          // The plan here is to make the new wedge design as a
          // "staircase".  That is, an extrusion whose face looks like
          // a series of steps.

          // First, get the extrusion HALF length
          double exHL = pbarWindow.width()*CLHEP::mm/2.0;

          // Now the thickness of each step
          double stepThck = pbarWindow.stripThickness()*CLHEP::mm;

          // The number of steps
          int nSteps = pbarWindow.nStrips();

          // The height of the strips, which is like the length of the
          // steps if you look at it as steps.
          std::vector<double> stepLength = pbarWindow.heights();

          //
          // the thickness of the steps, if variable
          std::vector<double> stepThickStrip;
          if (pbarAbsTS3Version==4) {
            stepThickStrip = pbarWindow.stripThicknesses();
          }

          // Sanity check.  There is no way we should get to this point
          // and have a stepLength vector with the wrong number of
          // steps, so just check...

          if ( (unsigned int) nSteps != stepLength.size() ) {
            throw cet::exception("GEOM")<<
              " The size of the PbarWedge stripHeight vector, "
                                        << stepLength.size() <<
              "\n Does not match the expected number of strips, "
                                        << nSteps << "\n" ;
          }

          if ( pbarAbsTS3Version == 4 && (unsigned int) nSteps != stepThickStrip.size() ) {
            throw cet::exception("GEOM")<<
              " The size of the PbarWedge stripThickness vector, "
                                        << stepThickStrip.size() <<
              "\n Does not match the expected number of strips, "
                                        << nSteps << "\n" ;
          }


          // Now we'll map out the vertices.  If you imagine our wedge
          // as a staircase leading up to a building, the origin of our
          // coordinate system will be the corner of the stairs on the bottom
          // and against the wall of the building.  The x-axis will point
          // up, against the wall, and the y-axis will run along the ground,
          // perpendicularly outward from the building.

          std::vector<G4TwoVector> stairOutline;
          stairOutline.reserve(2*nSteps + 2);  // # of vertices to describe

          // First point is (0,0)
          G4TwoVector tmpVertex(0,0);
          stairOutline.push_back(tmpVertex);
          double xCoord = 0.0;
          if (verbosityLevel > 2){
            G4cout << __func__ << "stepThck = " << stepThck << G4endl;
          }
          // Now loop over steps
          for ( int iStep = 0; iStep < nSteps; iStep++ ) {
            if ( verbosityLevel > 2) {
              G4cout << "istep, stair outline 1: " << iStep << " " << xCoord << " " << stepLength[iStep] << G4endl;
            }
            stairOutline.push_back(G4TwoVector(xCoord,-stepLength[iStep]));

            if (pbarAbsTS3Version == 3){
              xCoord += stepThck;
            } else  if (pbarAbsTS3Version == 4) {
              xCoord += stepThickStrip[iStep];
            }
            if ( verbosityLevel > 2) {
              G4cout << "istep, stair outline 2: " << iStep << " " << xCoord << " " << stepLength[iStep] << G4endl;
            }
            stairOutline.push_back(G4TwoVector(xCoord,-stepLength[iStep]));
          }

          G4TwoVector tmpVertex2(xCoord,0);
          stairOutline.push_back(tmpVertex2);

          G4ExtrudedSolid * stairCase = new G4ExtrudedSolid( "PbarAbsWedge_ES",
                                                             stairOutline,
                                                             exHL,
                                                             G4TwoVector(0,0),1.,
                                                             G4TwoVector(0,0),1. );


          AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
          G4RotationMatrix* pbarWedgeRot = reg.add(G4RotationMatrix());
          pbarWedgeRot->rotateY(90.0*CLHEP::degree);
          G4ThreeVector pbarWedgeTrans(0.0,pbarWedge_y1,pbarWedge_offsetZ);
          if (verbosityLevel > 2){
            G4cout << "pbarWedgeTrans = " << pbarWedgeTrans << G4endl;
          }
          pbarWedgeInfo.solid = stairCase;
          G4Material* wedgeMaterial  = findMaterialOrThrow(pbarWindow.wedgeMaterial()); //defaults to window material if not specified
          finishNesting(pbarWedgeInfo,
                        wedgeMaterial,
                        pbarWedgeRot,
                        pbarWedgeTrans,
                        parent.logical,
                        0,
                        G4Color::Yellow(),
                        "PbarAbs"
                        );
        } //end of if ( pbarWindow.shape == wedge )
    }  else {
      throw cet::exception("GEOM")<<
        " Incorrect pbar window geometry requested! \n " << " pbarWindow.shape() = " << pbarWindow.shape() << G4endl;

    } // end of else for pbarAbsTS3Version == ...


    // =============================================
    // ======== Now UPSTREAM pbar window ===========
    // =============================================

    // add a pbar window at the TS entrance
    TransportSolenoid const& ts = bl.getTS();
    GeomHandle<VirtualDetector> vdg;
    double vdHalfLength = vdg->getHalfLength()*CLHEP::mm;
    Mu2eG4Helper* _helper = &(*art::ServiceHandle<Mu2eG4Helper>() );

    bool is_pbarTS1In  = config.getBool("pbar.coll1In.build", true);
    bool is_pbarTS1Out = config.getBool("pbar.coll1Out.build", true);
    bool is_pbarTS31   = config.getBool("pbar.coll31In.build", false);

    // ************************************************************
    // Here we start building the upstream Pbar window.  Currently,
    // the model has it located upstream of the TS by 100 mm.
    // In the future we'd like to get this into the GeometryService
    // ***********************DNB**********************************

    if (is_pbarTS1In) {
      CollimatorTS1 const& coll1  = ts.getColl1() ;
      // ***
      // Pull in the existing parameters from config
      // ***
      string pbarTS1InMaterial   = config.getString("pbar.coll1In.material1Name");
      double pbarTS1InHalfLength = config.getDouble("pbar.coll1In.halfLength");
      double pbarTS1InROut       = config.getDouble("pbar.coll1In.rOut");
      double pbarTS1InRecordROut = config.getDouble("pbar.coll1In.rOutRecord");

      double pbarTS1InParams[5]  = { 0.0, pbarTS1InROut, pbarTS1InHalfLength, 0.0, CLHEP::twopi };  // coll1.rIn1()
      double pbarTS1InOffset = config.getDouble("pbar.coll1In.offset", 1.0);

      if ( verbosityLevel > 0 ) {
        G4cout << __func__ << " Pbar absorber at TS1 coll1 entrance halfLength : " << pbarTS1InHalfLength << G4endl;
        G4cout << __func__ << " Pbar absorber at TS1 coll1 entrance offset : " << pbarTS1InOffset << G4endl;
      }

      CLHEP::Hep3Vector pbarTS1InPos = coll1.getLocal();
      VolumeInfo motherVolume = _helper->locateVolInfo("TS1Vacuum");

      // This block determines the mother volume to use

      if (pbarTS1InOffset >= 0.0) {
        // use local when put in the TS1Vacuum
        pbarTS1InPos = coll1.getLocal();
        pbarTS1InPos.setZ( pbarTS1InPos.z() - coll1.halfLength() + 2.*vdHalfLength + pbarTS1InHalfLength + pbarTS1InOffset);
        motherVolume = _helper->locateVolInfo("TS1Vacuum");
      }
      else { // pbarTS1InOffset < 0.0
        // use global when put in the HallAir

        Tube const & psVacuumParams  = GeomHandle<PSVacuum>()->vacuum();

        pbarTS1InPos = ts.getTSVacuum<StraightSection>(TransportSolenoid::TSRegion::TS1)->getGlobal();
        pbarTS1InPos.setZ( pbarTS1InPos.z() - ts.getTSVacuum<StraightSection>(TransportSolenoid::TSRegion::TS1)->getHalfLength() - pbarTS1InHalfLength + pbarTS1InOffset);
        motherVolume = _helper->locateVolInfo("PSVacuum");
        G4ThreeVector psVacuumOriginInMu2e = psVacuumParams.originInMu2e();
        pbarTS1InPos = pbarTS1InPos - psVacuumOriginInMu2e;

        if ( verbosityLevel > 0 ) {

          G4cout << __func__ << " straight section halflength " << ts.getTSVacuum<StraightSection>(TransportSolenoid::TSRegion::TS1)->getHalfLength() << G4endl;
          G4cout << __func__ << " pbarTS1InHalfLength " << pbarTS1InHalfLength << G4endl;
          G4cout << __func__ << " pbarTS1InOffset " << pbarTS1InOffset << G4endl;
          G4cout << __func__ << " pbarTS1InPos " << pbarTS1InPos << G4endl;
        }
      }

      // mother volume set in block above

      // Here we put in the actual window - true in all versions
      nestTubs( "PbarAbsTS1In",
                pbarTS1InParams,
                findMaterialOrThrow(pbarTS1InMaterial),
                0,
                pbarTS1InPos,
                motherVolume,
                0,
                G4Color::Yellow(),
                "PbarAbs"
                );

      // ***
      // Check the version.  If version 2 or above, get new parameters
      // for the Pbar window supports.  December 2015
      // ***
      int pbarTS1InVersion = config.getInt("pbar.coll1In.Version",1);
      if ( pbarTS1InVersion > 1 ) {
        // Support structure inner and outer radius, halflength, and material
        double pbarTS1InSupRIn = config.getDouble("pbar.coll1In.supportRIn");
        double pbarTS1InSupROut = config.getDouble("pbar.coll1In.supportROut");
        double pbarTS1InSupHLen = config.getDouble("pbar.coll1In.supportHLen");
        string pbarTS1InSupMaterial = config.getString("pbar.coll1In.supportMaterialName");
        // Offset from Pbar window in z, and vector built from it
        double pbarTS1InSupOffsetZ = config.getDouble("pbar.coll1In.supportOffsetZ");
        CLHEP::Hep3Vector pbarTS1InSupRelPos(0,0,pbarTS1InSupOffsetZ);

        // set up the tube parameters
        double pbarTS1InSuptParams[5]  = { pbarTS1InSupRIn,
          pbarTS1InSupROut,
          pbarTS1InSupHLen,
          0.0, CLHEP::twopi };


        // Now put in the support ring
        nestTubs( "PbarAbsTS1InSup",
                  pbarTS1InSuptParams,
                  findMaterialOrThrow(pbarTS1InSupMaterial),
                  0,
                  pbarTS1InPos+pbarTS1InSupRelPos,
                  motherVolume,
                  0,
                  G4Color::Yellow(),
                  "PbarAbs"
                  );

        // The frames between which the Pbar window is sandwiched are
        // the same basic size and shape as the support ring, but with
        // different material and offsets.  So build in the same way
        // First get the new material and offset.

        string pbarTS1InFrameMaterial = config.getString("pbar.coll1In.frameMaterialName");
        double pbarTS1InFrameOffsetZ = config.getDouble("pbar.coll1In.frameOffsetZ");
        CLHEP::Hep3Vector pbarTS1InFrameRelPos(0,0,pbarTS1InFrameOffsetZ);

        nestTubs( "PbarAbsTS1InFrameUp",
                  pbarTS1InSuptParams,
                  findMaterialOrThrow(pbarTS1InFrameMaterial),
                  0,
                  pbarTS1InPos-pbarTS1InFrameRelPos,
                  motherVolume,
                  0,
                  G4Color::Yellow(),
                  "PbarAbs"
                  );


        nestTubs( "PbarAbsTS1InFrameDown",
                  pbarTS1InSuptParams,
                  findMaterialOrThrow(pbarTS1InFrameMaterial),
                  0,
                  pbarTS1InPos+pbarTS1InFrameRelPos,
                  motherVolume,
                  0,
                  G4Color::Yellow(),
                  "PbarAbs"
                  );

        // ***
        // Add the tabs used to hold the support ring in place.  Tabs are
        // located at 4, 8, and 12-o'clock
        // ***
        // Treat tab as a box.  Get half dimensions and material
        std::vector<double> pbarTS1InTabDims;
        config.getVectorDouble("pbar.coll1In.tabDims",pbarTS1InTabDims,3);
        string pbarTS1InTabMaterial = config.getString("pbar.coll1In.tabMaterialName");
        double pbarTS1InTabOffsetZ = config.getDouble("pbar.coll1In.tabOffsetZ");
        double pbarTS1InTabOffsetRad = config.getDouble("pbar.coll1In.tabOffsetR");

        CLHEP::Hep3Vector pbarTS1InTab1RelPos(0,
                                              pbarTS1InTabOffsetRad,
                                              pbarTS1InSupOffsetZ
                                              + pbarTS1InTabOffsetZ);
        // first tab at 12 o'clock
        nestBox( "PbarAbsTS1InSupTab1",
                 pbarTS1InTabDims,
                 findMaterialOrThrow(pbarTS1InTabMaterial),
                 0,
                 pbarTS1InPos+pbarTS1InTab1RelPos,
                 motherVolume,
                 0,
                 G4Color::Yellow(),
                 "PbarAbs"
                 );

        // next tab at 4 o'clock (looking downstream)
        CLHEP::Hep3Vector pbarTS1InTab2RelPos(-pbarTS1InTabOffsetRad*std::sin(120*CLHEP::degree),
                                              pbarTS1InTabOffsetRad*std::cos(120*CLHEP::degree),
                                              pbarTS1InSupOffsetZ+pbarTS1InTabOffsetZ);

        CLHEP::HepRotation* rotaTab2 = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
        rotaTab2->rotateZ(-120*CLHEP::degree);

        nestBox( "PbarAbsTS1InSupTab2",
                 pbarTS1InTabDims,
                 findMaterialOrThrow(pbarTS1InTabMaterial),
                 rotaTab2,
                 pbarTS1InPos+pbarTS1InTab2RelPos,
                 motherVolume,
                 0,
                 G4Color::Yellow(),
                 "PbarAbs"
                 );

        // next tab at 8 o'clock (looking downstream)
        CLHEP::Hep3Vector pbarTS1InTab3RelPos(pbarTS1InTabOffsetRad*std::sin(120*CLHEP::degree),
                                              pbarTS1InTabOffsetRad*std::cos(120*CLHEP::degree),
                                              pbarTS1InSupOffsetZ+pbarTS1InTabOffsetZ);

        CLHEP::HepRotation* rotaTab3 = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
        rotaTab3->rotateZ(120*CLHEP::degree);

        nestBox( "PbarAbsTS1InSupTab3",
                 pbarTS1InTabDims,
                 findMaterialOrThrow(pbarTS1InTabMaterial),
                 rotaTab3,
                 pbarTS1InPos+pbarTS1InTab3RelPos,
                 motherVolume,
                 0,
                 G4Color::Yellow(),
                 "PbarAbs"
                 );

        // ****
        // Now build the pegs used for handling.
        // Treat as simple polycones.
        // The pegs are located at 2, 6, and 10 o'clock around the z-axis
        // ****

        // Get shape and material parameters
        std::vector<double> pegIR;
        config.getVectorDouble("pbar.coll1In.pegInnerRadii",pegIR,7);
        std::vector<double> pegOR;
        config.getVectorDouble("pbar.coll1In.pegOuterRadii",pegOR,7);
        std::vector<double> pegZ;
        config.getVectorDouble("pbar.coll1In.pegZPlanes",pegZ,7);
        std::string pegMaterial = config.getString("pbar.coll1In.pegMaterialName");
        double pegOffsetZ = config.getDouble("pbar.coll1In.pegOffsetZ");
        double pegRad = config.getDouble("pbar.coll1In.pegRadialPosition");

        //This one at 2 o'clock looking down the z axis
        CLHEP::Hep3Vector peg1Pos(-pegRad*std::sin(60*CLHEP::degree),
                                  pegRad*std::cos(60*CLHEP::degree),
                                  -pegOffsetZ);

        Polycone pegCone( pegZ, pegIR, pegOR,
                          pbarTS1InPos - peg1Pos,
                          pegMaterial);

        VolumeInfo pegInfo = nestPolycone( "PbarAbsTS1InPeg1",
                                           pegCone.getPolyconsParams(),
                                           findMaterialOrThrow(pegMaterial),
                                           0,
                                           pbarTS1InPos+peg1Pos,
                                           motherVolume,
                                           0,
                                           G4Colour::Yellow(),
                                           "PbarAbs"
                                           );

        //This one at 6 o'clock looking down the z axis
        CLHEP::Hep3Vector peg2Pos(0, -pegRad, -pegOffsetZ);
        VolumeInfo peg2Info = nestPolycone( "PbarAbsTS1InPeg2",
                                            pegCone.getPolyconsParams(),
                                            findMaterialOrThrow(pegMaterial),
                                            0,
                                            pbarTS1InPos+peg2Pos,
                                            motherVolume,
                                            0,
                                            G4Colour::Yellow(),
                                            "PbarAbs"
                                            );


        //This one at 10 o'clock looking down the z axis
        CLHEP::Hep3Vector peg3Pos(pegRad*std::sin(60*CLHEP::degree),
                                  pegRad*std::cos(60*CLHEP::degree),
                                  -pegOffsetZ);

        VolumeInfo peg3Info = nestPolycone( "PbarAbsTS1InPeg3",
                                            pegCone.getPolyconsParams(),
                                            findMaterialOrThrow(pegMaterial),
                                            0,
                                            pbarTS1InPos+peg3Pos,
                                            motherVolume,
                                            0,
                                            G4Colour::Yellow(),
                                            "PbarAbs"
                                            );


      } // end of building support structure




      double pbarTS1InRecordParams[5]  = { 0.0, pbarTS1InRecordROut, vdHalfLength, 0.0, CLHEP::twopi };
      CLHEP::Hep3Vector pbarTS1InRecordPos = pbarTS1InPos;
      pbarTS1InRecordPos.setZ(pbarTS1InPos.z() - pbarTS1InHalfLength - 2*vdHalfLength - pbarTS1InRecordParams[2]);

      if ( verbosityLevel > 0 ) {
        G4cout << __func__ << "pbarTS1InRecordParams " << pbarTS1InRecordParams[1] << "  " << pbarTS1InRecordParams[2] << G4endl;
        G4cout << __func__ << "pbarTS1InRecordPos " << pbarTS1InRecordPos << G4endl;
      }

      nestTubs( "PbarAbsTS1InRecord",
                pbarTS1InRecordParams,
                findMaterialOrThrow(ts.upstreamVacuumMaterial()),
                0,
                pbarTS1InRecordPos,
                motherVolume,
                0,
                G4Color::Yellow(),
                "PbarAbs"
                );
    }

    if (is_pbarTS1Out) {

      // Get VDs
      // GeomHandle<VirtualDetector> vdg;
      // double vdHalfLength = vdg->getHalfLength()*CLHEP::mm;

      CollimatorTS1 const& coll1  = ts.getColl1() ;

      string pbarTS1OutMaterial   = config.getString("pbar.coll1Out.material1Name");
      // double pbarTS1OutHalfLength = config.getDouble("pbar.coll1Out.halfLength", 0.05);
      // double pbarTS1OutHalfLength = coll1.collarHalfLength()-2.*vdHalfLength;
      double pbarTS1OutHalfLength = coll1.collarHalfLength();
      double pbarTS1OutrIn        = coll1.collarrIn();
      double pbarTS1OutphiBegin   = coll1.collarphiBegin();
      double pbarTS1OutphiDelta   = coll1.collarphiDelta();
      double pbarTS1OutParams[5]  = { pbarTS1OutrIn, coll1.rIn1(), pbarTS1OutHalfLength,
        pbarTS1OutphiBegin*CLHEP::degree, pbarTS1OutphiDelta*CLHEP::degree };

      double pbarTS1OutPosz       = coll1.collarZ();

      if ( verbosityLevel > 0 ) {
        G4cout << __func__ << " Pbar absorber at TS1 coll1 near exit halfLength : " << pbarTS1OutHalfLength << " rIn " << pbarTS1OutrIn
               << " pbarTS1OutPosz " << pbarTS1OutPosz << " phiBegin " << pbarTS1OutphiBegin << " dPhi " << pbarTS1OutphiDelta << G4endl;
      }

      CLHEP::Hep3Vector pbarTS1OutPos = coll1.getLocal(); // the absorber is placed on the axis of coll1, shifted in z, see below
      // CLHEP::Hep3Vector TS1VacuumPos = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS1,TransportSolenoid::TSRadialPart::OUT)->getGlobal()-_hallOriginInMu2e;
      // pbarTS1OutPos.setZ( pbarTS1Outz - TS1VacuumPos.z() );

      if ( verbosityLevel > 1 ) {
        // printout related to the old code below
        G4cout << __func__ << " pbarTS1OutPos :                      " << pbarTS1OutPos << G4endl;
        G4cout << __func__ << " pbarTS1OutPosz-(-4044)) :            " << pbarTS1OutPosz-(-4044) << G4endl;
        G4cout << __func__ << " coll1.halfLength() :                 " << coll1.halfLength() << G4endl;
      }

      // pbarTS1OutPos.setZ( pbarTS1OutPos.z() + (pbarTS1OutPosz-(-4044)) - coll1.halfLength() );
      pbarTS1OutPos.setZ( coll1.collarZ() - ((_helper->locateVolInfo("TS1Vacuum")).centerInMu2e()).z() );

      if ( verbosityLevel > 0 ) {
        G4cout << __func__ << " PbarAbsTS1Out position in TS1Vacuum: " << pbarTS1OutPos << G4endl;
        double zpos = pbarTS1OutPos.z()+((_helper->locateVolInfo("TS1Vacuum")).centerInMu2e()).z();
        G4cout << __func__ << " PbarAbsTS1Out position in mu2e:      " << pbarTS1OutPos+(_helper->locateVolInfo("TS1Vacuum")).centerInMu2e() << G4endl;
        G4cout << __func__ << " PbarAbsTS1Out Extent:               [" << zpos - pbarTS1OutHalfLength << ","
               << zpos + pbarTS1OutHalfLength << "]" << G4endl;
        G4cout << __func__ << " PbarAbsTS1Out HalfLength :           " << pbarTS1OutHalfLength << G4endl;
      }

      nestTubs( "PbarAbsTS1Out",
                pbarTS1OutParams,
                findMaterialOrThrow(pbarTS1OutMaterial),
                0,
                pbarTS1OutPos,
                _helper->locateVolInfo("TS1Vacuum"),
                0,
                G4Color::Yellow(),
                "PbarAbs"
                );
    }

    if (is_pbarTS31) {
      CollimatorTS3 const& coll31 = ts.getColl31();

      double pbarTS31HalfLength = config.getDouble("pbar.coll31In.halfLength", 0.05);
      double pbarTS31Params[5]  = { 0.0, coll31.rOut(), pbarTS31HalfLength, 0.0, CLHEP::twopi };
      double pbarTS31Offset = config.getDouble("pbar.coll31In.offset", 1.0);

      CLHEP::Hep3Vector pbarTS31Pos = coll31.getLocal();
      pbarTS31Pos.setZ( pbarTS31Pos.z() - coll31.halfLength() - pbarTS31HalfLength - pbarTS31Offset);

      AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

      G4RotationMatrix* coll31Rot = reg.add(G4RotationMatrix());
      coll31Rot->rotateZ(coll31.rotationAngle()*CLHEP::degree);

      nestTubs( "PbarAbsTS31",
                pbarTS31Params,
                pbarMaterial,
                coll31Rot,
                pbarTS31Pos,
                _helper->locateVolInfo("TS3Vacuum"),
                0,
                G4Color::Yellow(),
                "PbarAbs"
                );
    }

  } // end Mu2eWorld::constructPbarWindow()

  //Helper function to add thermal shielding
  void addThermalShield ( TransportSolenoid const& ts, VolumeInfo const& useAsParent,
                          SimpleConfig const& config,
                          TransportSolenoid::TSRegion::enum_type TSRegion,
                          G4Material* thermalShieldMLIMaterial, G4Material* thermalShieldMidMaterial,
                          double centerWallThickness /* > 0 if exists in middle of region (TS3) */) {
    //for convenience
    typedef TransportSolenoid::TSRegion tsr;
    typedef TransportSolenoid::TSRadialPart tsrad;
    std::string region = ""; //store the number of the TS region
    bool straightSection = true; //torus or straight region
    switch(TSRegion) {
    case tsr::TS1 : region = "1"; break;
    case tsr::TS2 : region = "2"; straightSection = false; break;
    case tsr::TS3 : region = "3"; break;
    case tsr::TS4 : region = "4"; straightSection = false; break;
    case tsr::TS5 : region = "5"; break;
    default : return;
    }

    double smallGap = 1.e-3; //small buffer to prevent overlaps
    if(straightSection) {
      auto strsec = ts.getTSThermalShield<StraightSection>(TSRegion, tsrad::IN);
      //end plate for the thermal shielding
      double  endPlateRIn        = 0.;
      double  endPlateROut       = 0.;
      double  endPlateHalfLength = 0.;
      int side = (TSRegion == tsr::TS1) ? -1. : 1.; //which side of TS to add end plate to if adding it
      if(TSRegion == tsr::TS1 || TSRegion == tsr::TS5) { //only the TS1 and TS5 currently have plates
        endPlateRIn          = config.getDouble("ts.ts" + region + ".thermalshield.endplate.rIn");
        endPlateROut         = config.getDouble("ts.ts" + region + ".thermalshield.endplate.rOut");
        endPlateHalfLength   = config.getDouble("ts.ts" + region + ".thermalshield.endplate.halfLength");
        double endPlateZMid1 = config.getDouble("ts.ts" + region + ".thermalshield.endplate.mid.z1");
        double endPlateZMid2 = config.getDouble("ts.ts" + region + ".thermalshield.endplate.mid.z2");
        //calculate thickness of each blanket from total thickness - metal layer thickness
        if(endPlateZMid1 > endPlateZMid2 || endPlateZMid2 > 2.* endPlateHalfLength)
          throw cet::exception("GEOM")
            << ("TS thermal shielding end plate in TS" + region + " has inconsistent Z positions!").c_str()
            << " Should have 0 < z mid1 < z mid2 < 2*end plate half thickness."
            << G4endl;

        double midHalfLength = (endPlateZMid2 - endPlateZMid1)/2.; //z1 - z2
        double mli1HalfLength = endPlateZMid1/2.; //0 - z1
        double mli2HalfLength = endPlateHalfLength - endPlateZMid2/2.; //z2 - end
        double dzOriginToEdge = strsec->getHalfLength() - 2.*endPlateHalfLength; //center vacuum to edge of thermal shielding tube

        CLHEP::Hep3Vector mli1Origin(strsec->getGlobal()-useAsParent.centerInMu2e()
                                     + CLHEP::Hep3Vector(0.,0.,
                                                         side*(dzOriginToEdge + mli1HalfLength)));
        CLHEP::Hep3Vector midOrigin(mli1Origin + CLHEP::Hep3Vector(0.,0.,side*(mli1HalfLength+midHalfLength)));
        CLHEP::Hep3Vector mli2Origin(midOrigin + CLHEP::Hep3Vector(0.,0.,side*(mli2HalfLength+midHalfLength)));

        if(TSRegion == tsr::TS1){

          AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
          // Get the beam pipe
          G4Tubs* beamPassTub = nullptr;
          CLHEP::HepRotation* turn = nullptr;
          CLHEP::Hep3Vector place;
          getBeamPipe(config, reg, beamPassTub, turn, place);

          //
          std::string tssName1  = "TS"+region+"ThermalShieldEndPlateMLI1";
          G4VSolid *tssolid1= new G4Tubs( "ThermalShieldEndPlate",
                                          endPlateRIn ,endPlateROut,mli1HalfLength,0.,2*M_PI);
          CLHEP::Hep3Vector pos( mli1Origin+useAsParent.centerInWorld);
          G4SubtractionSolid* aSolid1 =
            new G4SubtractionSolid ( "ThermalShieldEndPlateMLI1",
                                     tssolid1,
                                     beamPassTub,
                                     turn,
                                     place-pos);
          VolumeInfo tss1(tssName1,
                          mli1Origin,
                          useAsParent.centerInWorld);
          tss1.solid = aSolid1;
          tss1.solid->SetName(tss1.name);

          finishNesting(tss1,
                        thermalShieldMLIMaterial,
                        0,
                        tss1.centerInParent,
                        useAsParent.logical,
                        0,
                        true,
                        G4Color::Blue(),
                        true,
                        true,
                        true,
                        false
                        );
          //
          std::string tssName2  = "TS"+region+"ThermalShieldEndPlateMLI2";
          G4VSolid *tssolid2= new G4Tubs( "ThermalShieldEndPlate",
                                          endPlateRIn ,endPlateROut,mli2HalfLength,0.,2*M_PI);
          pos= mli2Origin+useAsParent.centerInWorld;
          G4SubtractionSolid* aSolid2 =
            new G4SubtractionSolid ( "ThermalShieldEndPlateMLI2",
                                     tssolid2,
                                     beamPassTub,
                                     turn,
                                     place-pos);
          VolumeInfo tss2(tssName2,
                          mli2Origin,
                          useAsParent.centerInWorld);
          tss2.solid = aSolid2;
          tss2.solid->SetName(tss2.name);

          finishNesting(tss2,
                        thermalShieldMLIMaterial,
                        0,
                        tss2.centerInParent,
                        useAsParent.logical,
                        0,
                        true,
                        G4Color::Blue(),
                        true,
                        true,
                        true,
                        false
                        );
          //
          std::string tssNameMid  = "TS"+region+"ThermalShieldEndPlateMid";
          G4VSolid *tssolidMid= new G4Tubs( "ThermalShieldEndPlate",
                                            endPlateRIn ,endPlateROut,midHalfLength,0.,2*M_PI);

          pos= midOrigin+useAsParent.centerInWorld;
          G4SubtractionSolid* aSolidMid =
            new G4SubtractionSolid ( "ThermalShieldEndPlateMid",
                                     tssolidMid,
                                     beamPassTub,
                                     turn,
                                     place-pos);
          VolumeInfo tssMid(tssNameMid,
                            midOrigin,
                            useAsParent.centerInWorld);
          tssMid.solid = aSolidMid;
          tssMid.solid->SetName(tssMid.name);

          finishNesting(tssMid,
                        thermalShieldMidMaterial,
                        0,
                        tssMid.centerInParent,
                        useAsParent.logical,
                        0,
                        true,
                        G4Color::Red(),
                        true,
                        true,
                        true,
                        false
                        );

        } else {

          nestTubs( "TS"+region+"ThermalShieldEndPlateMLI1",
                    TubsParams( endPlateRIn ,
                                endPlateROut,
                                mli1HalfLength),
                    thermalShieldMLIMaterial,
                    0,
                    mli1Origin,
                    useAsParent,
                    0,
                    G4Color::Red(),
                    "TSCryo"
                    );
          nestTubs( "TS"+region+"ThermalShieldEndPlateMLI2",
                    TubsParams( endPlateRIn ,
                                endPlateROut,
                                mli2HalfLength),
                    thermalShieldMLIMaterial,
                    0,
                    mli2Origin,
                    useAsParent,
                    0,
                    G4Color::Red(),
                    "TSCryo"
                    );
          nestTubs( "TS"+region+"ThermalShieldEndPlateMid",
                    TubsParams( endPlateRIn ,
                                endPlateROut,
                                midHalfLength),
                    thermalShieldMidMaterial,
                    0,
                    midOrigin,
                    useAsParent,
                    0,
                    G4Color::Red(),
                    "TSCryo"
                    );
        }
      }
      if(centerWallThickness < 0.) { //no wall in the center of the region (TS1 and TS5)
        std::vector<double> innerInRadii = {strsec->rIn(),
          config.getDouble("ts.ts" + region + "in.thermalshield.mid.rIn")+smallGap,
          config.getDouble("ts.ts" + region + "in.thermalshield.mid.rOut")+smallGap};
        std::vector<double> innerOutRadii = {innerInRadii[1]-smallGap,
          innerInRadii[2]-smallGap,
          strsec->rOut()};
        CLHEP::Hep3Vector innerOrigin(strsec->getGlobal()-useAsParent.centerInMu2e()
                                      - CLHEP::Hep3Vector(0.,0.,side*endPlateHalfLength));
        double innerHalfLength = strsec->getHalfLength() - endPlateHalfLength;
        addThermalShieldStraightSection(config, useAsParent, innerInRadii, innerOutRadii, innerHalfLength,
                                        thermalShieldMLIMaterial, thermalShieldMidMaterial,
                                        innerOrigin, "TS"+region+"Inner");

        strsec = ts.getTSThermalShield<StraightSection>(TSRegion,tsrad::OUT);
        std::vector<double> outerInRadii = {strsec->rIn(),
          config.getDouble("ts.ts" + region + "out.thermalshield.mid.rIn")+smallGap,
          config.getDouble("ts.ts" + region + "out.thermalshield.mid.rOut")+smallGap};
        std::vector<double> outerOutRadii = {outerInRadii[1]-smallGap,
          outerInRadii[2]-smallGap,
          strsec->rOut()};
        CLHEP::Hep3Vector outerOrigin(strsec->getGlobal()-useAsParent.centerInMu2e()
                                      - CLHEP::Hep3Vector(0.,0.,side*endPlateHalfLength));
        double outerHalfLength = strsec->getHalfLength() - endPlateHalfLength;
        addThermalShieldStraightSection(config, useAsParent, outerInRadii, outerOutRadii, outerHalfLength,
                                        thermalShieldMLIMaterial, thermalShieldMidMaterial,
                                        outerOrigin, "TS"+region+"Outer");
      } //end no center wall
      else { //wall in the center of the region (TS3)
        //split each volume into an upstream and downstream end
        double tsrHalfLength = strsec->getHalfLength();
        double tsrudHalfLength = tsrHalfLength/2. - centerWallThickness/4.; //half of region minus half of the wall half length
        std::vector<double> innerInRadii = {strsec->rIn(),
          config.getDouble("ts.ts" + region + "in.thermalshield.mid.rIn")+smallGap,
          config.getDouble("ts.ts" + region + "in.thermalshield.mid.rOut")+smallGap};
        std::vector<double> innerOutRadii = {innerInRadii[1]-smallGap,
          innerInRadii[2]-smallGap,
          strsec->rOut()};
        CLHEP::Hep3Vector inneruOrigin(strsec->getGlobal()-useAsParent.centerInMu2e()
                                       - CLHEP::Hep3Vector(0.,0.,tsrHalfLength/2.+centerWallThickness/4.));
        CLHEP::Hep3Vector innerdOrigin(strsec->getGlobal()-useAsParent.centerInMu2e()
                                       + CLHEP::Hep3Vector(0.,0.,tsrHalfLength/2.+centerWallThickness/4.));

        addThermalShieldStraightSection(config, useAsParent, innerInRadii, innerOutRadii, tsrudHalfLength,
                                        thermalShieldMLIMaterial, thermalShieldMidMaterial,
                                        inneruOrigin, "TS"+region+"uInner");
        addThermalShieldStraightSection(config, useAsParent, innerInRadii, innerOutRadii, tsrudHalfLength,
                                        thermalShieldMLIMaterial, thermalShieldMidMaterial,
                                        innerdOrigin, "TS"+region+"dInner");

        strsec = ts.getTSThermalShield<StraightSection>(TSRegion,tsrad::OUT);
        tsrHalfLength = strsec->getHalfLength();
        tsrudHalfLength = tsrHalfLength/2. - centerWallThickness/4.; //half of region minus half of the wall half length
        std::vector<double> outerInRadii = {strsec->rIn(),
          config.getDouble("ts.ts" + region + "out.thermalshield.mid.rIn")+smallGap,
          config.getDouble("ts.ts" + region + "out.thermalshield.mid.rOut")+smallGap};
        std::vector<double> outerOutRadii = {outerInRadii[1]-smallGap,
          outerInRadii[2]-smallGap,
          strsec->rOut()};
        CLHEP::Hep3Vector outeruOrigin(strsec->getGlobal()-useAsParent.centerInMu2e()
                                       - CLHEP::Hep3Vector(0.,0.,tsrHalfLength/2.+centerWallThickness/4.));
        CLHEP::Hep3Vector outerdOrigin(strsec->getGlobal()-useAsParent.centerInMu2e()
                                       + CLHEP::Hep3Vector(0.,0.,tsrHalfLength/2.+centerWallThickness/4.));

        addThermalShieldStraightSection(config, useAsParent, outerInRadii, outerOutRadii, tsrudHalfLength,
                                        thermalShieldMLIMaterial, thermalShieldMidMaterial,
                                        outeruOrigin, "TS"+region+"uOuter");
        addThermalShieldStraightSection(config, useAsParent, outerInRadii, outerOutRadii, tsrudHalfLength,
                                        thermalShieldMLIMaterial, thermalShieldMidMaterial,
                                        outerdOrigin, "TS"+region+"dOuter");


      }
    } //end if straight section
    else { //torus section
      auto torsec = ts.getTSThermalShield<TorusSection>(TSRegion,tsrad::IN );
      std::vector<double> torusInParams = {torsec->torusRadius(), torsec->phiStart(), torsec->deltaPhi()};
      std::vector<double> innerInRadii = {torsec->rIn(),
        config.getDouble("ts.ts" + region + "in.thermalshield.mid.rIn"),
        config.getDouble("ts.ts" + region + "in.thermalshield.mid.rOut")};
      std::vector<double> innerOutRadii = {innerInRadii[1],
        innerInRadii[2],
        torsec->rOut()};
      addThermalShieldTorusSection(config, useAsParent, innerInRadii, innerOutRadii, torusInParams,
                                   thermalShieldMLIMaterial, thermalShieldMidMaterial,
                                   torsec->getGlobal()-useAsParent.centerInMu2e(), "TS"+region+"Inner");

      torsec = ts.getTSThermalShield<TorusSection>(TSRegion,tsrad::OUT );
      std::vector<double> torusOutParams = {torsec->torusRadius(), torsec->phiStart(), torsec->deltaPhi()};
      std::vector<double> outerInRadii = {torsec->rIn(),
        config.getDouble("ts.ts" + region + "out.thermalshield.mid.rIn"),
        config.getDouble("ts.ts" + region + "out.thermalshield.mid.rOut")};
      std::vector<double> outerOutRadii = {outerInRadii[1],
        outerInRadii[2],
        torsec->rOut()};
      addThermalShieldTorusSection(config, useAsParent, outerInRadii, outerOutRadii, torusOutParams,
                                   thermalShieldMLIMaterial, thermalShieldMidMaterial,
                                   torsec->getGlobal()-useAsParent.centerInMu2e(), "TS"+region+"Outer");
    }
  } //end addThermalShield

  //add an individual straight section of thermal shielding
  void addThermalShieldStraightSection( SimpleConfig const & config, VolumeInfo const& useAsParent, std::vector<double> innerRadii,
                                        std::vector<double> outerRadii, double halfLength,
                                        G4Material* thermalShieldMLIMaterial, G4Material* thermalShieldMidMaterial,
                                        CLHEP::Hep3Vector const& origin, std::string name) {
    unsigned nRadii = 3; //number expected
    if(innerRadii.size() != nRadii || outerRadii.size() != nRadii)
      throw cet::exception("GEOM")
        << "Incorrect amount of radii given for TS thermal shielding!"
        << G4endl;

    std::vector<std::string> names = {"ThermalShieldMLI1", "ThermalShieldMid", "ThermalShieldMLI2"};
    std::vector<G4Material*> materials = {thermalShieldMLIMaterial, thermalShieldMidMaterial, thermalShieldMLIMaterial};
    if(names.size() != nRadii)
      throw cet::exception("GEOM")
        << "Vector of shielding names not the correct size in constructTS::"
        << __func__ << "!"
        << G4endl;
    if(materials.size() != nRadii)
      throw cet::exception("GEOM")
        << "Vector of materials not the correct size in constructTS::"
        << __func__ << "!"
        << G4endl;
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    // Get the beam pipe
    G4Tubs* beamPassTub = nullptr;
    CLHEP::HepRotation* turn = nullptr;
    CLHEP::Hep3Vector place;
    if (name == "TS1Outer") getBeamPipe(config, reg, beamPassTub, turn, place);

    for(unsigned index = 0; index < nRadii; ++index) {
      if (name == "TS1Outer"){

        G4VSolid *tssolid= new G4Tubs( name+names[index], innerRadii[index], outerRadii[index], halfLength,0.,2*M_PI);

        CLHEP::Hep3Vector pos( origin+useAsParent.centerInWorld);
        G4SubtractionSolid* aSolid =
          new G4SubtractionSolid ( "TS1OuterFull",
                                   tssolid,
                                   beamPassTub,
                                   turn,
                                   place-pos);

        VolumeInfo tss(name+names[index],
                       origin,
                       useAsParent.centerInWorld);
        tss.solid = aSolid;
        tss.solid->SetName(tss.name);

        finishNesting(tss,
                      materials[index],
                      0,
                      tss.centerInParent,
                      useAsParent.logical,
                      0,
                      true,
                      G4Color::Blue(),
                      true,
                      true,
                      true,
                      false
                      );
      } else {
        nestTubs( name+names[index],
                  TubsParams( innerRadii[index],
                              outerRadii[index],
                              halfLength ),
                  materials[index],
                  0,
                  origin,
                  useAsParent,
                  0,
                  G4Color::Red(),
                  "TSCryo"
                  );
      }
    }
  } //end add thermal shield straight section

  void addThermalShieldTorusSection( SimpleConfig const & config, VolumeInfo const& useAsParent, std::vector<double> innerRadii,
                                     std::vector<double> outerRadii, std::vector<double> torusParams,
                                     G4Material* thermalShieldMLIMaterial, G4Material* thermalShieldMidMaterial,
                                     CLHEP::Hep3Vector const& origin, std::string name) {
    unsigned nRadii = 3; //number expected
    if(innerRadii.size() != nRadii || outerRadii.size() != nRadii)
      throw cet::exception("GEOM")
        << "Incorrect amount of radii given for TS thermal shielding!"
        << G4endl;
    if(torusParams.size() != 3)
      throw cet::exception("GEOM")
        << "Incorrect amount of torus parameters given for TS thermal shielding (expect {rTorus, phi0, delta phi})!"
        << G4endl;

    std::vector<std::string> names = {"ThermalShieldMLI1", "ThermalShieldMid", "ThermalShieldMLI2"};
    std::vector<G4Material*> materials = {thermalShieldMLIMaterial, thermalShieldMidMaterial, thermalShieldMLIMaterial};
    if(names.size() != nRadii)
      throw cet::exception("GEOM")
        << "Vector of shielding names not the correct size in constructTS::"
        << __func__ << "!"
        << G4endl;
    if(materials.size() != nRadii)
      throw cet::exception("GEOM")
        << "Vector of materials not the correct size in constructTS::"
        << __func__ << "!"
        << G4endl;

    //
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    // Get the beam pipe
    G4Tubs* beamPassTub = nullptr;
    CLHEP::HepRotation* turn = nullptr;
    CLHEP::Hep3Vector place;
    //
    for(unsigned index = 0; index < nRadii; ++index) {
      std::array<double,5> params {
        {   innerRadii[index], outerRadii[index],
            torusParams[0], torusParams[1],
            torusParams[2] }
      };
      if ( name == "TS2Outer" ){
        getBeamPipe(config, reg, beamPassTub, turn, place);

        G4VSolid *tssolid= new G4Torus( "TS2OuterFull",innerRadii[index], outerRadii[index], torusParams[0], torusParams[1], torusParams[2]);
        CLHEP::Hep3Vector pos( origin+useAsParent.centerInWorld);
        G4SubtractionSolid* aSolid =
          new G4SubtractionSolid ( "TS2OuterFull",
                                   tssolid,
                                   beamPassTub,
                                   turn,
                                   place-pos);
        VolumeInfo tss(name+names[index],
                       origin,
                       useAsParent.centerInWorld);
        tss.solid = aSolid;
        tss.solid->SetName(tss.name);
        finishNesting(tss,
                      materials[index],
                      0,
                      tss.centerInParent,
                      useAsParent.logical,
                      0,
                      true,
                      G4Color::Blue(),
                      true,
                      true,
                      true,
                      false
                      );
      }
      else{
        nestTorus(name+names[index],
                  params,
                  materials[index],
                  0,
                  origin,
                  useAsParent,
                  0,
                  G4Color::Red(),
                  "TSCryo"
                  );
      }
    }

  }//end add thermal shield torus section

  // Function to get beam pipe dimensions, orientation and position
  void getBeamPipe(SimpleConfig const& config, AntiLeakRegistry& reg, G4Tubs* &beamPassTub, CLHEP::HepRotation * &turn,CLHEP::Hep3Vector &place){
    GeomHandle<PSShield> hrs;
    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "PSShield", "PSShield" );
    Tube const & pssInletParams = hrs->beamInlet();

    // Define the beam inlet asa tube: use 2mm as safety margin in Rmax and 4100 mm as additional half length in Z to cover both PS and TS
    beamPassTub = new G4Tubs( "beampipePassthrough",
                              0.0,
                              pssInletParams.getTubsParams().data()[1]+2*CLHEP::mm,
                              pssInletParams.getTubsParams().data()[2]+4100.0*CLHEP::mm,
                              0.0*CLHEP::degree,
                              360.0*CLHEP::degree );

    turn = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
    turn->rotateY(hrs->getBeamAngleY());
    turn->rotateX(hrs->getBeamAngleX());

    place = hrs->getBeamInletCenter();

  }
} //end namespace mu2e
