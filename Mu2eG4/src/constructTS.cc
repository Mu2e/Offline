//
// Free function to create Transport Solenoid
//
// $Id: constructTS.cc,v 1.33 2014/09/19 19:15:10 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/09/19 19:15:10 $
//
// Original author KLG based on Mu2eWorld constructTS
//
// Notes:
// Construct the TS.  Parent volume is the air inside of the hall.

// Mu2e includes.
#include "Mu2eG4/inc/constructTS.hh"
#include "G4Helper/inc/VolumeInfo.hh"

// C++ includes
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Framework includes
#include "cetlib_except/exception.h"

// Mu2e includes.
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/Collimator_TS1.hh"
#include "BeamlineGeom/inc/PbarWindow.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "GeomPrimitives/inc/PolyconsParams.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestTorus.hh"
#include "Mu2eG4/inc/nestCons.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "GeomPrimitives/inc/Polycone.hh"
#include "ProductionSolenoidGeom/inc/PSVacuum.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"

using namespace std;

namespace mu2e {

  void constructTS( VolumeInfo const & parent,
                    SimpleConfig const & c
                    ) {

    GeomHandle<Beamline> bl;
    G4Helper* _helper = &(*art::ServiceHandle<G4Helper>() );

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

    G4Helper * const _helper = &(*art::ServiceHandle<G4Helper>());

    const int  verbosityLevel      = config.getInt ("ts.cryo.verbosityLevel", 0     );
    const bool polyLiningUp        = config.getBool("ts.polyliner.Up.build"   , false );
    const bool polyLiningDown      = config.getBool("ts.polyliner.Down.build" , false );

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "TSCryo", "ts.cryo"      );
    geomOptions->loadEntry( config, "TSPoly", "ts.polyliner" );
    geomOptions->loadEntry( config, "TSCA"  , "ts.cas"       );
    // For how all pieces are made from one of two types of material,
    // vacuum or average coils + cryostat material.
    G4Material* downstreamVacuumMaterial  = findMaterialOrThrow(ts->downstreamVacuumMaterial());
    G4Material* upstreamVacuumMaterial    = findMaterialOrThrow(ts->upstreamVacuumMaterial());
    G4Material* cryoMaterial    = findMaterialOrThrow(ts->material());

    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

    // Build upstream end wall of TS1
    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS1,TransportSolenoid::TSRadialPart::IN );

    CLHEP::Hep3Vector pos( strsec->getGlobal().x(), 
                           strsec->getGlobal().y(), 
                           strsec->getGlobal().z()-strsec->getHalfLength()+ts->endWallU1_halfLength() );

    std::string tssName  = "TS1UpstreamEndwall";
    nestTubs( tssName,
              TubsParams( ts->endWallU1_rIn(),
                          ts->endWallU1_rOut(),
                          ts->endWallU1_halfLength() ),
              cryoMaterial,
              0,
              pos-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
	      "TSCryo"
              );

    verbosityLevel &&
      G4cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << G4endl;
    if ( verbosityLevel ) {
      G4cout << __func__ << " Upstream TS1 endwall at: " << pos << G4endl;
    }

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
      G4cout << __func__ << " TS1(out) OffsetInMu2e  : " << strsec->getGlobal()   << G4endl;
      G4cout << __func__ << " TS1(out) rotation      : " << *(strsec->getRotation()) << G4endl;
    }
    
    // Build downstream partial end wall of TS1
    CLHEP::Hep3Vector pos2( strsec->getGlobal().x(), 
                            strsec->getGlobal().y(), 
                            strsec->getGlobal().z()+strsec->getHalfLength()+ts->endWallU2_halfLength() );
  
    tssName = "TS1DownstreamEndwall";
    nestTubs( tssName,
              TubsParams( ts->endWallU2_rIn(),
                          ts->endWallU2_rOut(),
                          ts->endWallU2_halfLength() ),
              cryoMaterial,
              0,
              pos2-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Red(),
	      "TSCryo"
              );
              
    verbosityLevel &&
      G4cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << G4endl;

    if ( verbosityLevel ) {
      G4cout << __func__ << " Downstream TS1 endwall at: " << pos2 << G4endl;
      G4cout << __func__ << " Downstream TS1 extent   [: " << pos2.z()-ts->endWallU2_halfLength() 
           << "," << pos2.z() +ts->endWallU2_halfLength() << "]" << G4endl;
    }

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
      VolumeInfo ts2Vacuum = art::ServiceHandle<G4Helper>()->locateVolInfo("TS2Vacuum");
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
    std::array<double,5> ts2Cryo2Params { { torsec->rIn(), torsec->rOut(), torsec->torusRadius(), torsec->phiStart(), torsec->deltaPhi() } };

    tssName = "TS2OuterCryoShell";
    double ts2InsVacROut = torsec->rIn();
    nestTorus(tssName,
              ts2Cryo2Params,
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


    // Build TS3
    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS3,TransportSolenoid::TSRadialPart::IN );
    nestTubs( "TS3Vacuum",
              TubsParams( 0., ts->innerRadius(), strsec->getHalfLength() ),
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
      VolumeInfo ts4Vacuum = art::ServiceHandle<G4Helper>()->locateVolInfo("TS4Vacuum");
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




    // Build TS5.
    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS5,TransportSolenoid::TSRadialPart::IN );    
    CLHEP::Hep3Vector globalVac5Position = CLHEP::Hep3Vector(strsec->getGlobal().x(),
                                                             strsec->getGlobal().y(),
                                                             strsec->getGlobal().z() );

    nestTubs( "TS5Vacuum",
              TubsParams( 0., ts->innerRadius(), strsec->getHalfLength() ),
              downstreamVacuumMaterial,
              strsec->getRotation(),
              globalVac5Position-_hallOriginInMu2e,
              parent,
              0,
              G4Color::Green(),
	      "TSCryo"
              );

    double ts5InsVacRIn = strsec->rOut();
    tssName = "TS5InnerCryoShell";
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
      CLHEP::HepRotation* ringRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
      double ringRotTheta = thetasRing[iRing]*CLHEP::degree;
      ringRotat->rotateY(ringRotTheta);
      CLHEP::HepRotation* noRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);

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
      leftName << "leftSideRing" << iRing;

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
      centerName << "centerRing" << iRing;

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
      rightName << "rightSideRing" << iRing;

      double rx = 0.0;
      double ry = 0.0;
      double rz = lr/2.0 + trs/2.0; 

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

    } // End of building rings

    // Build outer portion of support for TS3 antiproton stopping window
    // Will model for now as solid tubes

    int pbarAbsTS3Version = config.getInt("pbar.version",1);
      if ( pbarAbsTS3Version > 1 ) {

	std::ostringstream PabsSupOutName;
	PabsSupOutName << "PabsTS3SupOut";
	CLHEP::HepRotation* pasubRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
	pasubRotat->rotateY(90.0*CLHEP::degree);

	nestTubs( PabsSupOutName.str(),
		  TubsParams( rirs, rors, trs ),
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

	G4Helper* _helper = &(*art::ServiceHandle<G4Helper>() );
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
      G4Helper* _helper = &(*art::ServiceHandle<G4Helper>() );

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
      } else if ( its==tsCAReg_enum::TS1 
                  || its==tsCAReg_enum::TS3u 
                  || its==tsCAReg_enum::TS3d 
                  || its==tsCAReg_enum::TS5
                  ) {

	if ( its == tsCAReg_enum::TS1 ) {
	  useAsParent = _helper->locateVolInfo( "TS1CryoInsVac" );
	} else if ( its == tsCAReg_enum::TS3u || its == tsCAReg_enum::TS3d ) {
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

	G4Helper* _helper = &(*art::ServiceHandle<G4Helper>() );
	VolumeInfo useAsParent = _helper->locateVolInfo( ivName.str() );

	CLHEP::Hep3Vector V = coil.getGlobal() - useAsParent.centerInMu2e();


	CLHEP::HepRotation* rot = useAsParent.physical->GetRotation();
	V = (*rot) * V;
	CLHEP::HepRotation invrot = rot->inverse();

	if ( iTS == TransportSolenoid::TSRegion::TS2 || 
	     iTS == TransportSolenoid::TSRegion::TS4 ) {

	  CLHEP::HepRotation* twist = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
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

    TSSection const * ts1in = ts.getTSCryo(TransportSolenoid::TSRegion::TS1, TransportSolenoid::TSRadialPart::IN);
    TSSection const * ts3in = ts.getTSCryo(TransportSolenoid::TSRegion::TS3, TransportSolenoid::TSRadialPart::IN);
    TSSection const * ts5in = ts.getTSCryo(TransportSolenoid::TSRegion::TS5, TransportSolenoid::TSRadialPart::IN);

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

    G4Helper* _helper = &(*art::ServiceHandle<G4Helper>() );

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
      G4cout << __func__ << " ts innerRadius      : " << ts.innerRadius() << G4endl;
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
    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

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
    TubsParams coll31OutRecordParam ( 0,  ts.innerRadius(), vdHalfLength );
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


    TubsParams coll5Param2 ( coll52.rIn(),  coll52.rOut(), coll52.halfLength()-2.*vdHalfLength);
    nestTubs( "Coll52",
              coll5Param2,
              findMaterialOrThrow( coll52.material() ),
              0,
              coll52.getLocal(),
              _helper->locateVolInfo("TS5Vacuum"),
              0,
              G4Color::Blue(),
	      "TSColl"
              );
    
    TubsParams coll5Param3 ( coll53.rIn(),  coll53.rOut(), coll53.halfLength()-2.*vdHalfLength);
    
    // Build Coll53 flange outside the TS5Vacuum
    VolumeInfo useAsParent = _helper->locateVolInfo("TS5CryoInsVac");
    //    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();
    G4ThreeVector _useParentOriginInMu2e = useAsParent.centerInMu2e();
    StraightSection   const * strsec (nullptr);
    strsec = ts.getTSVacuum<StraightSection>(TransportSolenoid::TSRegion::TS5);
    nestTubs( "Coll53",
              coll5Param3,
              findMaterialOrThrow( coll53.material() ),
	      0,
              strsec->getGlobal()-_useParentOriginInMu2e+coll53.getLocal(),
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

        AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
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
      
	  AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
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

      CLHEP::HepRotation * supportRot = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
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
	// nestTubs( "PbarAbs",
	// 	  pbarParams,
	// 	  pbarMaterial,
	// 	  0,
	// 	  pbarWindow.getLocal(),
	// 	  parent,
	// 	  0,
	// 	  G4Color::Yellow()
	// 	  );
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

	  AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
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

      CLHEP::HepRotation * supportRot = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);

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


	  AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
	  G4RotationMatrix* pbarWedgeRot = reg.add(G4RotationMatrix());
	  pbarWedgeRot->rotateY(90.0*CLHEP::degree);
	  G4ThreeVector pbarWedgeTrans(0.0,pbarWedge_y1,pbarWedge_offsetZ);
	  if (verbosityLevel > 2){
	    G4cout << "pbarWedgeTrans = " << pbarWedgeTrans << G4endl;
	  }
	  pbarWedgeInfo.solid = stairCase;
	  
	  finishNesting(pbarWedgeInfo,
			pbarMaterial,
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
    G4Helper* _helper = &(*art::ServiceHandle<G4Helper>() );

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

	CLHEP::HepRotation* rotaTab2 = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
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

	CLHEP::HepRotation* rotaTab3 = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
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

      AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

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

}

