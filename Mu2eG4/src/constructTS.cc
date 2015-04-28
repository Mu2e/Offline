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

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Framework includes
#include "cetlib/exception.h"

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
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestTorus.hh"
#include "Mu2eG4/inc/nestCons.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeomPrimitives/inc/Tube.hh"
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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;
    if ( verbosityLevel ) {
      cout << __func__ << " Upstream TS1 endwall at: " << pos << endl;
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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;
    
    if ( verbosityLevel > 0) {
      cout << __func__ << " TS1(in)  OffsetInMu2e  : " << strsec->getGlobal()   << endl;
      cout << __func__ << " TS1(in)  Extent        :[ " << strsec->getGlobal().z() - strsec->getHalfLength() <<","  
           << strsec->getGlobal().z() + strsec->getHalfLength() << "]" << endl;
      cout << __func__ << " TS1(in)  rotation      : " << strsec->getRotation() << endl;
     }

    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS1,TransportSolenoid::TSRadialPart::OUT );
    tssName =  "TS1OuterCryoShell";
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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;

    if ( verbosityLevel > 0) {
      cout << __func__ << " TS1(out) OffsetInMu2e  : " << strsec->getGlobal()   << endl;
      cout << __func__ << " TS1(out) rotation      : " << strsec->getRotation() << endl;
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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;

    if ( verbosityLevel ) {
      cout << __func__ << " Downstream TS1 endwall at: " << pos2 << endl;
      cout << __func__ << " Downstream TS1 extent   [: " << pos2.z()-ts->endWallU2_halfLength() 
           << "," << pos2.z() +ts->endWallU2_halfLength() << "]" << endl;
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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;

    torsec = ts->getTSCryo<TorusSection>(TransportSolenoid::TSRegion::TS2,TransportSolenoid::TSRadialPart::OUT );
    std::array<double,5> ts2Cryo2Params { { torsec->rIn(), torsec->rOut(), torsec->torusRadius(), torsec->phiStart(), torsec->deltaPhi() } };

    tssName = "TS2OuterCryoShell";
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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;

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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;

    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS3,TransportSolenoid::TSRadialPart::OUT );

    tssName = "TS3OuterCryoShell";
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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;

    if ( verbosityLevel > 0) {
      cout << __func__ << " TS3  OffsetInMu2e : " << strsec->getGlobal()   << endl;
      cout << __func__ << " TS3  rotation     : " << strsec->getRotation() << endl;
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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;

    torsec = ts->getTSCryo<TorusSection>(TransportSolenoid::TSRegion::TS4,TransportSolenoid::TSRadialPart::OUT );
    std::array<double,5> ts4Cryo2Params { { torsec->rIn(), torsec->rOut(), torsec->torusRadius(), torsec->phiStart(), torsec->deltaPhi() } };
    tssName = "TS4OuterCryoShell";
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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;

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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;

    strsec = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS5,TransportSolenoid::TSRadialPart::OUT );
    tssName = "TS5OuterCryoShell";
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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;
   
    if ( verbosityLevel > 0) {
      cout << __func__ << " TS5  OffsetInMu2e : " << strsec->getGlobal()   << endl;
      cout << __func__ << " TS5  rotation     : " << strsec->getRotation() << endl;
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
    
    for ( unsigned int iRing = 0; iRing < xr.size(); iRing++ ) {
      std::ostringstream leftName;
      leftName << "leftSideRing" << iRing;
      CLHEP::HepRotation* ringRotat = new CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY);
      double ringRotTheta = 0.0;
      if ( iRing == 1 || iRing == 4 ) ringRotTheta = 45.0*CLHEP::degree;
      if ( iRing == 2 || iRing == 3 ) ringRotTheta = 90.0*CLHEP::degree;
      ringRotat->rotateY(ringRotTheta);
      double lx = xr[iRing] + lr*sin(ringRotTheta)/2.0 + trs*sin(ringRotTheta)/2.0;
      double ly = yr[iRing];
      double lz = zr[iRing] - lr*cos(ringRotTheta)/2.0 - trs*cos(ringRotTheta)/2.0; 
      nestTubs( leftName.str(),
		TubsParams( rirs, rors, trs/2.0 ),
		ringMaterial,
                ringRotat,
		CLHEP::Hep3Vector(lx,ly,lz)-_hallOriginInMu2e,
		parent,
		0,
		G4Color::Blue(),
		"TSCryo"
		);

      std::ostringstream centerName;
      centerName << "centerRing" << iRing;

      nestTubs( centerName.str(),
		TubsParams( rir, ror, lr/2.0 ),
		ringMaterial,
                ringRotat,
		CLHEP::Hep3Vector(xr[iRing],yr[iRing],zr[iRing])-_hallOriginInMu2e,
		parent,
		0,
		G4Color::Blue(),
		"TSCryo"
		);


      std::ostringstream rightName;
      rightName << "rightSideRing" << iRing;

      double rx = xr[iRing] - lr*sin(ringRotTheta)/2.0 - trs*sin(ringRotTheta)/2.0;
      double ry = yr[iRing];
      double rz = zr[iRing] + lr*cos(ringRotTheta)/2.0 + trs*cos(ringRotTheta)/2.0; 

      nestTubs( rightName.str(),
		TubsParams( rirs, rors, trs/2.0 ),
		ringMaterial,
                ringRotat,
		CLHEP::Hep3Vector(rx,ry,rz)-_hallOriginInMu2e,
		parent,
		0,
		G4Color::Blue(),
		"TSCryo"
		);

    }

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
      std::cout << __func__ << " " << tssName << " Mass in kg: " 
                << _helper->locateVolInfo(tssName).logical->GetMass()/CLHEP::kg 
                << std::endl;

    if ( verbosityLevel ) {
      cout << __func__ << " Downstream TS5 endwall at: " << pos3 << endl;
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
      verbosityLevel && std::cout << __func__ << " constructing " << caName << std::endl;

      if ( its==tsCAReg_enum::TS2 || its==tsCAReg_enum::TS4 ) {
        
        // those sections are toruses
        // fixme, make the base class TSSection more general to be used here instead the two
        caTorsec = ts.getTSCA<TorusSection>(its);
        nestTorus(caName,
                  caTorsec->getParameters(),
                  findMaterialOrThrow(caTorsec->getMaterial()),
                  caTorsec->getRotation(),
                  caTorsec->getGlobal()-parentCenterInMu2e,
                  parent,
                  0,
                  G4Color::Yellow(),
                  "TSCA"
                  );

      } else if ( its==tsCAReg_enum::TS1 
                  || its==tsCAReg_enum::TS3u 
                  || its==tsCAReg_enum::TS3d 
                  || its==tsCAReg_enum::TS5
                  ) {

        caConsec = ts.getTSCA<ConeSection>(its);
        nestCons(caName,
                 caConsec->getParameters(),
                 findMaterialOrThrow(caConsec->getMaterial()),
                 caConsec->getRotation(),
                 caConsec->getGlobal()-parentCenterInMu2e,
                 parent,
                 0,
                 G4Color::Yellow(),
                 "TSCA"
                 );

      } else {

        caStrsec = ts.getTSCA<StraightSection>(its);
        nestTubs(caName,
                 TubsParams(caStrsec->rIn(),
                            caStrsec->rOut(),
                            caStrsec->getHalfLength()),
                 findMaterialOrThrow(caStrsec->getMaterial()),
                 caStrsec->getRotation(),
                 caStrsec->getGlobal()-parentCenterInMu2e,
                 parent,
                 0,
                 G4Color::Cyan(),
                 "TSCA"
                 );

      }

      G4Helper* _helper = &(*art::ServiceHandle<G4Helper>());
      verbosityLevel 
        && std::cout << __func__ << " " << caName << " Mass in kg: " 
                     << _helper->locateVolInfo(caName).logical->GetMass()/CLHEP::kg 
                     << std::endl;
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

        nestTubs( coilname.str(),
                  TubsParams( coil.rIn(), coil.rOut(), coil.halfLength() ),
                  coilMaterial,
                  coil.getRotation(),
                  coil.getGlobal()-parent.centerInMu2e(),
                  parent,
                  0,
                  G4Color::Green(),
		  "TSCoils"
                  );

        if ( verbosityLevel > 0 ) {
          cout << __func__ << " " << coilname.str() << " placed at: " << coil.getGlobal() << endl;
          cout << __func__ << "            rotation: " << -coil.getRotation()->getTheta()/CLHEP::degree << endl;
          cout << __func__ << "              params: " << coil.rIn() << " , " << coil.rOut() << " , " << 2*coil.halfLength() << endl;
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

    CollimatorTS1 const& coll1  = ts.getColl1() ;
    CollimatorTS3 const& coll31 = ts.getColl31();
    CollimatorTS3 const& coll32 = ts.getColl32();
    CollimatorTS5 const& coll5  = ts.getColl5() ;

    // Get VDs
    GeomHandle<VirtualDetector> vdg;
    double vdHalfLength = vdg->getHalfLength()*CLHEP::mm;

    // Place collimator 1 (concentric cone which can be a cylinder when r1==r2
    // and a cylinder placed in TS1Vacuum)

    // the cone (which can be a tube/cylinder) inside the outer tube/cylinder
    double coll1Param1[7] = { coll1.rIn1(), coll1.rIn3(),
                              coll1.rIn2(), coll1.rIn3(),
                              coll1.halfLength() - 2.*vdHalfLength,
                              0.0, CLHEP::twopi };

    G4Helper* _helper = &(*art::ServiceHandle<G4Helper>() );

    CLHEP::Hep3Vector parentPos    = _helper->locateVolInfo("TS1Vacuum").centerInWorld;

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

    TubsParams coll1Param2 ( coll1.rIn3(),  ts.innerRadius(), coll1.halfLength()-2.*vdHalfLength);
 
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
    
    if ( verbosityLevel > 0) {
      cout << __func__ << " TS1  OffsetInMu2e    : " << ts1in->getGlobal()       << endl;
      cout << __func__ << " Coll1 local offset   : " << ts.getColl1().getLocal() << endl;
      cout << __func__ << " TS1  Rotation        : " << ts1in->getRotation()     << endl;
    }

    // Place collimator 3

    // Collimator 3 has peculiar shape, described in doc_db 853.
    // Construct this shape using boolean functions on solids

    // First, construct hole; make it slightly longer that any collimator
    double hDz = coll31.halfLength();
    if( hDz<coll32.halfLength() ) hDz=coll32.halfLength();
    // Hole is the intersection of box and circles
    G4Box* coll3_hole_box = new G4Box("coll3_hole_box",
                                      coll31.holeRadius()+5.0,coll31.holeHalfHeight(),hDz+1.0);
    G4Tubs* coll3_hole_circle = new G4Tubs("coll3_hole_circle",
                                           0.0,coll31.holeRadius(),hDz+1.0,
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

    G4Tubs* coll31_mother = new G4Tubs("Coll31_mother",
                                       0, ts.innerRadius(), coll31.halfLength()-2.*vdHalfLength,
                                       0.0, CLHEP::twopi );

    G4Tubs* coll32_mother = new G4Tubs("Coll32_mother",
                                       0, ts.innerRadius(), coll32.halfLength()-2.*vdHalfLength,
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

    TubsParams coll32InRecordParam ( 0,  ts.innerRadius(), vdHalfLength );
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
      cout << __func__ << " TS3  OffsetInMu2e   : " << ts3in->getGlobal() << endl;
      cout << __func__ << " Coll31 local offest : " << coll31.getLocal() << endl;
      cout << __func__ << " Coll32 local offset : " << coll32.getLocal() << endl;
      cout << __func__ << " TS3  Rotation       : " << ts3in->getRotation() << endl;
    }

    // Place collimator 5

    if ( verbosityLevel > 0) {
      cout << __func__ << " TS5  OffsetInMu2e  : " << ts5in->getGlobal()   << endl;
      cout << __func__ << " Coll5 local offset : " << coll5.getLocal()             << endl;
      cout << __func__ << " TS5  Rotation      : " << ts5in->getRotation() << endl;
    }

    CLHEP::Hep3Vector coll5OffsetInMu2e = ts5in->getGlobal() + 
      ( ( ts5in->getRotation() != 0x0 ) ?
        *(ts5in->getRotation()) * coll5.getLocal() : 
        coll5.getLocal() );

    if ( verbosityLevel > 0) {
      cout << __func__ << "  coll5OffsetInMu2e    : "    << coll5OffsetInMu2e << endl;
      cout << __func__ << "  Coll5 calc local offset : " << coll5OffsetInMu2e - ts5in->getGlobal() << endl;
    }

    // the most outer part (with Virtual Detectors on the outer surfaces of the Coll5)

    Tube coll5Param(coll5.material(),
                    coll5OffsetInMu2e,
                    coll5.rIn(),
                    coll5.rOut() - 2.*vdHalfLength,
                    coll5.halfLength() - 2.*vdHalfLength);

    VolumeInfo coll5Info = nestTubs( "Coll5",
                                     coll5Param.getTubsParams(),
                                     findMaterialOrThrow(coll5Param.materialName()),
                                     0,
                                     coll5.getLocal(),
                                     _helper->locateVolInfo("TS5Vacuum"),
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

    CollimatorTS5 const& coll5 = bl.getTS().getColl5();

    if( degraderR.size()!=degraderDZB.size() || degraderR.size()!=degraderDZT.size() ||
        degraderR.size()!=degraderPhi.size() ) {
      cout << __func__ << " Warning: MuonDegrader is not build - dimensions don't match." << endl;
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
          cout << __func__ << " Degrader constructed at: " << coll5.getLocal() << " wrt. TS5 " << endl;
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

    // Place Pbar absorber between Coll31 and Coll32
    // Pbar absorber is made of two pieces:
    //  -- vacuum wall, which covers the whole inner part of TS3
    //     it is controlled by pbar.* parameters
    //  -- wedge, which starts near center and extends upward
    //     it is controlled by pbarwedge.* parameters

    // -- vacuum wall
    int const verbosityLevel = config.getInt("pbar.verbosityLevel", 0);

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "PbarAbs", "pbar" );
    

    // Throw exception if pbarwedge.build is used
    if ( config.hasName("pbarwedge.build") )
      {
        throw cet::exception("GEOM")<<
          " Variable pbarwedge.build is now deprecated. \n" <<
          " To use pbar wedge specify: pbar.Type = \"wedge\" \n" ;
      }


    PbarWindow const & pbarWindow = bl.getTS().getPbarWindow();
    G4Material* pbarMaterial  = findMaterialOrThrow( pbarWindow.material() );
    if (verbosityLevel > 0) std::cout << "TS3 pbar windows HalfLength : " << pbarWindow.halfLength() << std::endl; 
        
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
        
    if( pbarWindow.shape() == "wedge" ) 
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
                      G4ThreeVector(0.,0.,pbarWedge_dz/2+pbarWindow.halfLength()),
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
    else 
      {
        throw cet::exception("GEOM")<<
          " Incorrect pbar window geometry requested! \n " ;
      }

    // add a pbar window at the TS entrance
    TransportSolenoid const& ts = bl.getTS();
    GeomHandle<VirtualDetector> vdg;
    double vdHalfLength = vdg->getHalfLength()*CLHEP::mm;
    G4Helper* _helper = &(*art::ServiceHandle<G4Helper>() );

    bool is_pbarTS1In  = config.getBool("pbar.coll1In.build", true);
    bool is_pbarTS1Out = config.getBool("pbar.coll1Out.build", true);
    bool is_pbarTS2    = config.getBool("pbar.TS2.build", true);
    bool is_pbarTS31   = config.getBool("pbar.coll31In.build", false);

    if (is_pbarTS1In) {
      CollimatorTS1 const& coll1  = ts.getColl1() ;

      string pbarTS1InMaterial   = config.getString("pbar.coll1In.material1Name");
      double pbarTS1InHalfLength = config.getDouble("pbar.coll1In.halfLength");
      double pbarTS1InROut       = config.getDouble("pbar.coll1In.rOut");
      double pbarTS1InRecordROut = config.getDouble("pbar.coll1In.rOutRecord");

      double pbarTS1InParams[5]  = { 0.0, pbarTS1InROut, pbarTS1InHalfLength, 0.0, CLHEP::twopi };  // coll1.rIn1()
      double pbarTS1InOffset = config.getDouble("pbar.coll1In.offset", 1.0);
      if ( verbosityLevel > 0 ) {
        cout << __func__ << "Pbar absorber at TS1 coll1 entrance halfLength : " << pbarTS1InHalfLength << std::endl;
        cout << __func__ << "Pbar absorber at TS1 coll1 entrance offset : " << pbarTS1InOffset << std::endl;
      }

      CLHEP::Hep3Vector pbarTS1InPos = coll1.getLocal();
      VolumeInfo motherVolume = _helper->locateVolInfo("TS1Vacuum");
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
          cout << __func__ << "coll1 halflength " << ts.getTSVacuum<StraightSection>(TransportSolenoid::TSRegion::TS1)->getHalfLength() << endl;
          cout << __func__ << "pbarTS1InHalfLength " << pbarTS1InHalfLength << endl;
          cout << __func__ << "pbarTS1InOffset " << pbarTS1InOffset << endl;
          cout << __func__ << "pbarTS1InPos " << pbarTS1InPos << endl;
        }
      }

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

      double pbarTS1InRecordParams[5]  = { 0.0, pbarTS1InRecordROut, vdHalfLength, 0.0, CLHEP::twopi };  
      CLHEP::Hep3Vector pbarTS1InRecordPos = pbarTS1InPos;
      pbarTS1InRecordPos.setZ(pbarTS1InPos.z() - pbarTS1InHalfLength - 2*vdHalfLength - pbarTS1InRecordParams[2]);

      if ( verbosityLevel > 0 ) {
        cout << __func__ << "pbarTS1InRecordParams " << pbarTS1InRecordParams[1] << "  " << pbarTS1InRecordParams[2] << endl;
        cout << __func__ << "pbarTS1InRecordPos " << pbarTS1InRecordPos << endl;
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
      CollimatorTS1 const& coll1  = ts.getColl1() ;

      string pbarTS1OutMaterial   = config.getString("pbar.coll1Out.material1Name");
      double pbarTS1OutHalfLength = config.getDouble("pbar.coll1Out.halfLength", 0.05);
      double pbarTS1OutrIn        = config.getDouble("pbar.coll1Out.rIn",        120.0);
      double pbarTS1OutphiBegin   = config.getDouble("pbar.coll1Out.phiBegin",   210.0);
      double pbarTS1OutphiDelta   = config.getDouble("pbar.coll1Out.phiDelta",   120.0);
      double pbarTS1OutParams[5]  = { pbarTS1OutrIn, coll1.rIn1(), pbarTS1OutHalfLength,
                                      pbarTS1OutphiBegin*CLHEP::degree, pbarTS1OutphiDelta*CLHEP::degree };
      double pbarTS1OutPosz       = config.getDouble("pbar.coll1Out.z", -3144.0);
      if ( verbosityLevel > 0 ) {
        std::cout << "Pbar absorber at TS1 coll1 near exit halfLength : " << pbarTS1OutHalfLength << " rIn " << pbarTS1OutrIn 
          << " pbarTS1OutPosz " << pbarTS1OutPosz << " phiBegin " << pbarTS1OutphiBegin << " dPhi " << pbarTS1OutphiDelta << std::endl;
      }

      CLHEP::Hep3Vector pbarTS1OutPos = coll1.getLocal();
//      CLHEP::Hep3Vector TS1VacuumPos = ts->getTSCryo<StraightSection>(TransportSolenoid::TSRegion::TS1,TransportSolenoid::TSRadialPart::OUT)->getGlobal()-_hallOriginInMu2e;
//      pbarTS1OutPos.setZ( pbarTS1Outz - TS1VacuumPos.z() );
      pbarTS1OutPos.setZ( pbarTS1OutPos.z() + (pbarTS1OutPosz-(-4044)) - coll1.halfLength() );

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

    if (is_pbarTS2) {
      double rTorus = config.getDouble("ts.rTorus", 2929.0);

      vector<double> pbarTS2Origin;
      config.getVectorDouble("pbar.TS2.origin", pbarTS2Origin, vector<double>() );
      CLHEP::Hep3Vector originTS2(pbarTS2Origin[0], pbarTS2Origin[1], pbarTS2Origin[2]); 

      double theta = config.getDouble("pbar.TS2.theta", 85.0);   
      CLHEP::Hep3Vector pbarTS2Pos(sin(theta*CLHEP::degree), -cos(theta*CLHEP::degree), 0.0);
      pbarTS2Pos *= rTorus;
      //cout << "pbarTS2 pos " << pbarTS2Pos << endl;

      double pbarTS2HalfLength = config.getDouble("pbar.TS2.halfLength", 0.5);
      double pbarTS2rIn        = config.getDouble("pbar.TS2.rIn",        100.0);
      double pbarTS2rOut       = config.getDouble("pbar.TS2.rOut",       240.0);
      double pbarTS2phiBegin   = config.getDouble("pbar.TS2.phiBegin",   210.0);
      double pbarTS2phiDelta   = config.getDouble("pbar.TS2.phiDelta",   120.0);
      double pbarTS2Params[5]  = { pbarTS2rIn, pbarTS2rOut, pbarTS2HalfLength,
                                   pbarTS2phiBegin*CLHEP::degree, pbarTS2phiDelta*CLHEP::degree };

      std::cout << "pbarTS2 params: halfLength " << pbarTS2HalfLength << " rIn " << pbarTS2rIn << " rOut " << pbarTS2rOut
           << " phiBegin " << pbarTS2phiBegin << " dPhi " << pbarTS2phiDelta << std::endl;
    
      AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
      G4RotationMatrix* pbarTS2Rot = reg.add(G4RotationMatrix());
      pbarTS2Rot->rotateX( -90.0*CLHEP::degree );
      pbarTS2Rot->rotateY( (90.0-theta)*CLHEP::degree );

      nestTubs( "PbarAbsTS2",
                pbarTS2Params,
                pbarMaterial,
                pbarTS2Rot,
                pbarTS2Pos,
                _helper->locateVolInfo("TS2Vacuum"),
                0,
                G4Color::Yellow(),
		"PbarAbs"
              );
    }

    if (is_pbarTS31) {
      CollimatorTS3 const& coll31 = ts.getColl31();

      double pbarTS31HalfLength = config.getDouble("pbar.coll31In.halfLength", 0.05);
      double pbarTS31Params[5]  = { 0.0, ts.innerRadius(), pbarTS31HalfLength, 0.0, CLHEP::twopi };
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

