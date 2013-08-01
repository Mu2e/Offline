//
// Free function to create Transport Solenoid
//
// $Id: constructTS.cc,v 1.23 2013/08/01 14:37:08 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/08/01 14:37:08 $
//
// Original author KLG based on Mu2eWorld constructTS
//
// Notes:
// Construct the TS.  Parent volume is the air inside of the hall.

// Mu2e includes.
#include "Mu2eG4/inc/constructTS.hh"
// #include "Mu2eG4/inc/ConstructTransportSolenoid.hh"
#include "G4Helper/inc/VolumeInfo.hh"

// C++ includes
#include <iostream>

// Framework includes
#include "cetlib/exception.h"

// Mu2e includes.
#include "BeamlineGeom/inc/Beamline.hh"
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

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"

using namespace std;

namespace mu2e {

  void constructTS( VolumeInfo const & parent,
                    SimpleConfig const & c
                    ) {

    GeomHandle<Beamline> bl;
    G4Helper* _helper = &(*art::ServiceHandle<G4Helper>() );

    constructCryostat   ( parent, c, *bl);
    constructCoils      ( parent, c, *bl);
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

    const int  verbosityLevel      = config.getInt("ts.cryo.verbosityLevel", 0);
    const bool visible             = config.getBool("ts.cryo.visible",true);
    const bool solid               = config.getBool("ts.cryo.solid",true);
    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV             = true;

    // For how all pieces are made from one of two types of material,
    // vacuum or average coils + cryostat material.
    G4Material* vacuumMaterial  = findMaterialOrThrow(ts->insideMaterial());
    G4Material* cryoMaterial    = findMaterialOrThrow(ts->material());

    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

    const double rTorus   = ts->torusRadius();
    const double rVac     = ts->innerRadius();

    // Build upstream end wall of TS1
    strsec = &ts->getTS1_in();

    CLHEP::Hep3Vector pos( strsec->getGlobal().x(), 
                           strsec->getGlobal().y(), 
                           strsec->getGlobal().z()-strsec->getHalfLength()-ts->endWallU1_halfLength() );
  
    nestTubs( "TS1UpstreamEndwall",
              TubsParams( ts->endWallU1_rIn(),
                          ts->endWallU1_rOut(),
                          ts->endWallU1_halfLength() ),
              cryoMaterial,
              0,
              pos-_hallOriginInMu2e,
              parent,
              0,
              visible,
              G4Color::Red(),
              solid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );
              
    if ( verbosityLevel ) {
      cout << __func__ << " Upstream TS1 endwall at: " << pos << endl;
    }

    // Build TS1.
    CLHEP::Hep3Vector globalPosition = CLHEP::Hep3Vector(strsec->getGlobal().x(),
                                                         strsec->getGlobal().y(),
                                                         strsec->getGlobal().z()-ts->endWallU1_halfLength() );
    VolumeInfo ts1VacInfo = nestTubs( "TS1Vacuum",
                                      TubsParams(0., rVac, 
                                                 strsec->getHalfLength() + ts->endWallU1_halfLength() ),
                                      vacuumMaterial,
                                      strsec->getRotation(),
                                      globalPosition-_hallOriginInMu2e,
                                      parent,
                                      0,
                                      visible,
                                      G4Color::Red(),
                                      solid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    VolumeInfo ts1CryoInfo1 = nestTubs( "TS1InnerCryoShell",
                                       TubsParams( strsec->rIn(),
                                                   strsec->rOut(),
                                                   strsec->getHalfLength() ),
                                       cryoMaterial,
                                       strsec->getRotation(),
                                       strsec->getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       visible,
                                       G4Color::Red(),
                                       solid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    if ( verbosityLevel > 0) {
      cout << __func__ << " TS1(in)  OffsetInMu2e  : " << strsec->getGlobal()   << endl;
      cout << __func__ << " TS1(in)  Extent        :[ " << strsec->getGlobal().z() - strsec->getHalfLength() + ts->endWallU1_halfLength() <<","  
           << strsec->getGlobal().z() + strsec->getHalfLength() - ts->endWallU1_halfLength() << "]" << endl;
      cout << __func__ << " TS1(in)  rotation      : " << strsec->getRotation() << endl;
    }

    strsec = &ts->getTS1_out();
    VolumeInfo ts1CryoInfo2 = nestTubs( "TS1OuterCryoShell",
                                       TubsParams( strsec->rIn(),
                                                   strsec->rOut(),
                                                   strsec->getHalfLength() ),
                                       cryoMaterial,
                                       strsec->getRotation(),
                                       strsec->getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       visible,
                                       G4Color::Red(),
                                       solid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    if ( verbosityLevel > 0) {
      cout << __func__ << " TS1(out) OffsetInMu2e  : " << strsec->getGlobal()   << endl;
      cout << __func__ << " TS1(out) rotation      : " << strsec->getRotation() << endl;
    }

    // Build downstream partial end wall of TS1
    CLHEP::Hep3Vector pos2( strsec->getGlobal().x(), 
                            strsec->getGlobal().y(), 
                            strsec->getGlobal().z()+strsec->getHalfLength()+ts->endWallU2_halfLength() );
  
    nestTubs( "TS1DownstreamEndwall",
              TubsParams( ts->endWallU2_rIn(),
                          ts->endWallU2_rOut(),
                          ts->endWallU2_halfLength() ),
              cryoMaterial,
              0,
              pos2-_hallOriginInMu2e,
              parent,
              0,
              visible,
              G4Color::Red(),
              solid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );
              
    if ( verbosityLevel ) {
      cout << __func__ << " Downstream TS1 endwall at: " << pos2 << endl;
      cout << __func__ << " Donwstream TS1 extent   [: " << pos2.z()-ts->endWallU1_halfLength() << "," << pos2.z() +ts->endWallU1_halfLength() << "]" << endl;
    }

    // Build TS2.
    torsec = &ts->getTS2_in();

    double ts2VacParams[5]   = { 0.0,   rVac, rTorus, 
                                 torsec->phiStart(), torsec->deltaPhi() };
    double ts2Cryo1Params[5] = { torsec->rIn(), torsec->rOut(), torsec->torusRadius(), 
                                 torsec->phiStart(), torsec->deltaPhi() };

    VolumeInfo ts2VacInfo = nestTorus("TS2Vacuum",
                                      ts2VacParams,
                                      vacuumMaterial,
                                      torsec->getRotation(),
                                      torsec->getGlobal()-_hallOriginInMu2e,
                                      parent,
                                      0,
                                      visible,
                                      G4Color::Green(),
                                      solid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    VolumeInfo ts2Cryo1Info = nestTorus("TS2InnerCryoShell",
                                        ts2Cryo1Params,
                                        cryoMaterial,
                                        torsec->getRotation(),
                                        torsec->getGlobal()-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        visible,
                                        G4Color::Red(),
                                        solid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );

    torsec = &ts->getTS2_out();
    double ts2Cryo2Params[5] = { torsec->rIn(), torsec->rOut(), torsec->torusRadius(), 
                                 torsec->phiStart(), torsec->deltaPhi() };

    VolumeInfo ts2Cryo2Info = nestTorus("TS2OuterCryoShell",
                                        ts2Cryo2Params,
                                        cryoMaterial,
                                        torsec->getRotation(),
                                        torsec->getGlobal()-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        visible,
                                        G4Color::Red(),
                                        solid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );

    // Build TS3.
    strsec = &ts->getTS3_in();
    VolumeInfo ts3VacInfo = nestTubs( "TS3Vacuum",
                                      TubsParams( 0., rVac, strsec->getHalfLength() ),
                                      vacuumMaterial,
                                      strsec->getRotation(),
                                      strsec->getGlobal()-_hallOriginInMu2e,
                                      parent,
                                      0,
                                      visible,
                                      G4Color::Green(),
                                      solid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    VolumeInfo ts3CryoInfo1 = nestTubs( "TS3InnerCryoShell",
                                       TubsParams( strsec->rIn(),
                                                   strsec->rOut(),
                                                   strsec->getHalfLength() ),
                                       cryoMaterial,
                                       strsec->getRotation(),
                                       strsec->getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       visible,
                                       G4Color::Red(),
                                       solid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    strsec = &ts->getTS3_out();
    VolumeInfo ts3CryoInfo2 = nestTubs( "TS3OuterCryoShell",
                                        TubsParams( strsec->rIn(),
                                                    strsec->rOut(),
                                                    strsec->getHalfLength() ),
                                        cryoMaterial,
                                        strsec->getRotation(),
                                        strsec->getGlobal()-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        visible,
                                        G4Color::Red(),
                                        solid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );

    if ( verbosityLevel > 0) {
      cout << __func__ << " TS3  OffsetInMu2e : " << strsec->getGlobal()   << endl;
      cout << __func__ << " TS3  rotatio      : " << strsec->getRotation() << endl;
    }

    // Build TS4.
    torsec = &ts->getTS4_in();
    double ts4VacParams[5]   = { 0.0,   rVac, rTorus, 
                                 torsec->phiStart(), torsec->deltaPhi() };
    double ts4Cryo1Params[5] = { torsec->rIn(), torsec->rOut(), torsec->torusRadius(), 
                                 torsec->phiStart(), torsec->deltaPhi() };

    VolumeInfo ts4VacInfo = nestTorus("TS4Vacuum",
                                      ts4VacParams,
                                      vacuumMaterial,
                                      torsec->getRotation(),
                                      torsec->getGlobal()-_hallOriginInMu2e,
                                      parent,
                                      0,
                                      visible,
                                      G4Color::Yellow(),
                                      solid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    VolumeInfo ts4CryoInfo1 = nestTorus("TS4InnerCryoShell",
                                       ts4Cryo1Params,
                                       cryoMaterial,
                                       torsec->getRotation(),
                                       torsec->getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       visible,
                                       G4Color::Red(),
                                       solid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    torsec = &ts->getTS4_out();
    double ts4Cryo2Params[5] = { torsec->rIn(), torsec->rOut(), torsec->torusRadius(), 
                                 torsec->phiStart(), torsec->deltaPhi() };

    VolumeInfo ts4CryoInfo2 = nestTorus("TS4OuterCryoShell",
                                       ts4Cryo2Params,
                                       cryoMaterial,
                                       torsec->getRotation(),
                                       torsec->getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       visible,
                                       G4Color::Red(),
                                       solid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );
    // Build TS5.
    strsec = &ts->getTS5_in();
    
    CLHEP::Hep3Vector globalVac5Position = CLHEP::Hep3Vector(strsec->getGlobal().x(),
                                                             strsec->getGlobal().y(),
                                                             strsec->getGlobal().z()+ts->endWallD_halfLength() );

    VolumeInfo ts5VacInfo = nestTubs( "TS5Vacuum",
                                      TubsParams( 0., rVac, 
                                                  strsec->getHalfLength()+ts->endWallD_halfLength() ),
                                      vacuumMaterial,
                                      strsec->getRotation(),
                                      globalVac5Position-_hallOriginInMu2e,
                                      parent,
                                      0,
                                      visible,
                                      G4Color::Green(),
                                      solid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    VolumeInfo ts5CryoInfo1 = nestTubs( "TS5InnerCryoShell",
                                       TubsParams( strsec->rIn(),
                                                   strsec->rOut(),
                                                   strsec->getHalfLength() ),
                                       cryoMaterial,
                                       strsec->getRotation(),
                                       strsec->getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       visible,
                                       G4Color::Red(),
                                       solid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    strsec = &ts->getTS5_out();
    VolumeInfo ts5CryoInfo2 = nestTubs( "TS5OuterCryoShell",
                                       TubsParams( strsec->rIn(),
                                                   strsec->rOut(),
                                                   strsec->getHalfLength() ),
                                       cryoMaterial,
                                       strsec->getRotation(),
                                       strsec->getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       visible,
                                       G4Color::Red(),
                                       solid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    
    if ( verbosityLevel > 0) {
      cout << __func__ << " TS5  OffsetInMu2e : " << strsec->getGlobal()   << endl;
      cout << __func__ << " TS5  rotatio      : " << strsec->getRotation() << endl;
    }

    // Build downstream end wall of TS5
    CLHEP::Hep3Vector pos3( strsec->getGlobal().x(), 
                            strsec->getGlobal().y(), 
                            strsec->getGlobal().z()+strsec->getHalfLength()+ts->endWallD_halfLength() );
  
    nestTubs( "TS5DownstreamEndwall",
              TubsParams( ts->endWallD_rIn(),
                          ts->endWallD_rOut(),
                          ts->endWallD_halfLength() ),
              cryoMaterial,
              0,
              pos3-_hallOriginInMu2e,
              parent,
              0,
              visible,
              G4Color::Red(),
              solid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    if ( verbosityLevel ) {
      cout << __func__ << " Downstream TS5 endwall at: " << pos3 << endl;
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
    const bool visible             = config.getBool("ts.coils.visible",true);
    const bool solid               = config.getBool("ts.coils.solid",false);
    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV             = true;

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
                  visible,
                  G4Color::Green(),
                  solid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

        if ( verbosityLevel > 0 ) {
          cout << __func__ << " " << coilname.str() << " placed at: " << coil.getGlobal() << endl;
          cout << __func__ << "            rotation: " << -coil.getRotation()->getTheta()/degree << endl;
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
    const bool visible             = config.getBool("ts.coll.visible",true);
    const bool solid               = config.getBool("ts.coll.solid",true);
    const int  verbosityLevel      = config.getInt("ts.coll.verbosityLevel", 0);
    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV             = true;

    // Get collimators
    TransportSolenoid const& ts = bl.getTS();

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

    nestCons( "Coll11",
              coll1Param1,
              findMaterialOrThrow( coll1.material1() ),
              0,
              coll1.getLocal(),
              _helper->locateVolInfo("TS1Vacuum"),
              0,
              visible,
              G4Color::Cyan(),
              solid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    TubsParams coll1Param2 ( coll1.rIn3(),  ts.innerRadius(), coll1.halfLength()-2.*vdHalfLength);
 
    nestTubs( "Coll12",
              coll1Param2,
              findMaterialOrThrow( coll1.material2() ),
              0,
              coll1.getLocal(),
              _helper->locateVolInfo("TS1Vacuum"),
              0,
              visible,
              G4Color::Blue(),
              solid,
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );
    
    if ( verbosityLevel > 0) {
      cout << __func__ << " TS1  OffsetInMu2e    : " << ts.getTS1_in().getGlobal()   << endl;
      cout << __func__ << " Coll5 local offset   : " << ts.getColl5().getLocal()     << endl;
      cout << __func__ << " TS1  Rotation        : " << ts.getTS1_in().getRotation() << endl;
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
                  visible,
                  G4Color::Gray(),
                  solid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);

    finishNesting(coll32Info,
                  findMaterialOrThrow( coll32.material() ),
                  coll32Rot,
                  coll32.getLocal(),
                  _helper->locateVolInfo("TS3Vacuum").logical,
                  0,
                  visible,
                  G4Color::Gray(),
                  solid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);

    if ( verbosityLevel > 0) {
      cout << __func__ << " TS3  OffsetInMu2e   : " << ts.getTS3_in().getGlobal() << endl;
      cout << __func__ << " Coll31 local offest : " << coll31.getLocal() << endl;
      cout << __func__ << " Coll32 local offset : " << coll32.getLocal() << endl;
      cout << __func__ << " TS3  Rotation       : " << ts.getTS3_in().getRotation() << endl;
    }

    // Place collimator 5

    if ( coll5.rOut() > ts.innerRadius() || 
         coll5.rMid2()  > coll5.rOut() - 2.*vdHalfLength || 
         coll5.rMid1()  >= coll5.rMid2() || 
         coll5.rIn() > coll5.rMid1() ||
         coll5.rIn() < 0.0 ) {

      throw cet::exception("GEOM")<< " constructTS: wrong coll5 radii: " 
                                  << "\n rVac             : " << ts.innerRadius()
                                  << "\n coll5OuterRadius : " << coll5.rOut()
                                  << "\n coll5MidRadius2  : " << coll5.rMid2()
                                  << "\n coll5MidRadius1  : " << coll5.rMid1()
                                  << "\n coll5InnerRadiu  : " << coll5.rIn()
                                  << "\n";
    }

    if ( coll5.halfLengthU() < 0.0 || coll5.halfLengthD() < 0.0 ||
         coll5.halfLengthU() + coll5.halfLengthD() > coll5.halfLength() - 2.*vdHalfLength) {

      throw cet::exception("GEOM")<< " constructTS: wrong coll5 longitudinal params " 
                                  << "\n coll5HalfLength   : " << coll5.halfLength()
                                  << "\n coll5HalfLengthU  : " << coll5.halfLengthU()
                                  << "\n coll5HalfLengthD  : " << coll5.halfLengthD()
                                  << "\n";
    }

    if ( verbosityLevel > 0) {
      cout << __func__ << " TS5  OffsetInMu2e  : " << ts.getTS5_in().getGlobal()   << endl;
      cout << __func__ << " Coll5 local offset : " << coll5.getLocal()             << endl;
      cout << __func__ << " TS5  Rotation      : " << ts.getTS5_in().getRotation() << endl;
    }

    CLHEP::Hep3Vector coll5OffsetInMu2e = ts.getTS5_in().getGlobal() + 
      ( ( ts.getTS5_in().getRotation() != 0x0 ) ?
        *(ts.getTS5_in().getRotation()) * coll5.getLocal() : 
        coll5.getLocal() );

    if ( verbosityLevel > 0) {
      cout << __func__ << "  coll5OffsetInMu2e    : "    << coll5OffsetInMu2e << endl;
      cout << __func__ << "  Coll5 calc local offset : " << coll5OffsetInMu2e - ts.getTS5_in().getGlobal() << endl;
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
                                     visible,
                                     G4Color::Blue(),
                                     solid,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

    // decide if we need to add absorber parts

    if ( coll5.halfLengthD() > 0.0 ) {
      
      // + downstream absorber
      
      CLHEP::Hep3Vector coll5AbsorberOffsetInColl5 = 
        CLHEP::Hep3Vector (0.0, 0.0, coll5.halfLength() - 2.*vdHalfLength - coll5.halfLengthD());

      CLHEP::Hep3Vector coll5AbsorberOffsetInMu2e = coll5OffsetInMu2e +
        ( ( ts.getTS5_in().getRotation() != 0x0 ) ?
          *(ts.getTS5_in().getRotation()) * coll5AbsorberOffsetInColl5 :
          coll5AbsorberOffsetInColl5 );

      Tube coll5DAbsParam(coll5.absMaterial(),
                          coll5AbsorberOffsetInMu2e,
                          coll5.rIn(),
                          coll5.rOut() - 2.*vdHalfLength,
                          coll5.halfLengthD());

      nestTubs( "coll5DAbsorber",
                coll5DAbsParam.getTubsParams(),
                findMaterialOrThrow(coll5DAbsParam.materialName()),
                0,
                coll5AbsorberOffsetInColl5,
                coll5Info,
                0,
                visible,
                G4Color::Gray(),
                solid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
    }

    if ( coll5.halfLengthU() > 0.0 ) {

      // + upstream absorber
    
      CLHEP::Hep3Vector coll5AbsorberOffsetInColl5 = 
        CLHEP::Hep3Vector (0.0, 0.0, - coll5.halfLength() + 2.*vdHalfLength + coll5.halfLengthU());

      CLHEP::Hep3Vector coll5AbsorberOffsetInMu2e = coll5OffsetInMu2e +
        ( ( ts.getTS5_in().getRotation() != 0x0 ) ?
          *(ts.getTS5_in().getRotation()) * coll5AbsorberOffsetInColl5 :
          coll5AbsorberOffsetInColl5 );

      Tube coll5AbsParam(coll5.absMaterial(),
                         coll5AbsorberOffsetInMu2e,
                         coll5.rIn(),
                         coll5.rOut() - 2.*vdHalfLength,
                         coll5.halfLengthU());

      nestTubs( "coll5UAbsorber",
                coll5AbsParam.getTubsParams(),
                findMaterialOrThrow(coll5AbsParam.materialName()),
                0,
                coll5AbsorberOffsetInColl5,
                coll5Info,
                0,
                visible,
                G4Color::Gray(),
                solid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
    }

    if ( coll5.rMid1() > coll5.rIn() ) {

      // + "inner radius" absorber
    
      CLHEP::Hep3Vector coll5AbsorberOffsetInColl5 =
        CLHEP::Hep3Vector (0.0, 0.0,
                           coll5.halfLengthU() - coll5.halfLengthD());

      CLHEP::Hep3Vector coll5AbsorberOffsetInMu2e = coll5OffsetInMu2e +
        ( ( ts.getTS5_in().getRotation() != 0x0 ) ?
          *(ts.getTS5_in().getRotation()) * coll5AbsorberOffsetInColl5 :
          coll5AbsorberOffsetInColl5 );

      Tube coll5AbsParam(coll5.absMaterial(),
                         coll5AbsorberOffsetInMu2e,
                         coll5.rIn(),
                         coll5.rMid1(),
                         coll5.halfLength() - 2.*vdHalfLength - 
                         coll5.halfLengthD() - coll5.halfLengthU());

      nestTubs( "coll5IAbsorber",
                coll5AbsParam.getTubsParams(),
                findMaterialOrThrow(coll5AbsParam.materialName()),
                0,
                coll5AbsorberOffsetInColl5,
                coll5Info,
                0,
                visible,
                G4Color::Gray(),
                solid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
    }

    if ( coll5.rMid2() < coll5.rOut() - 2.*vdHalfLength) {

      // + "outer radius" absorber
    
      CLHEP::Hep3Vector coll5AbsorberOffsetInColl5 =
        CLHEP::Hep3Vector (0.0, 0.0,
                           coll5.halfLengthU() - coll5.halfLengthD());

      CLHEP::Hep3Vector coll5AbsorberOffsetInMu2e = coll5OffsetInMu2e +
        ( ( ts.getTS5_in().getRotation() != 0x0 ) ?
          *(ts.getTS5_in().getRotation()) * coll5AbsorberOffsetInColl5 :
          coll5AbsorberOffsetInColl5 );

      Tube coll5AbsParam(coll5.absMaterial(),
                         coll5AbsorberOffsetInMu2e,
                         coll5.rMid2(),
                         coll5.rOut() - 2.*vdHalfLength,
                         coll5.halfLength() - 2.*vdHalfLength - 
                         coll5.halfLengthD() - coll5.halfLengthU());

      nestTubs( "coll5OAbsorber",
                coll5AbsParam.getTubsParams(),
                findMaterialOrThrow(coll5AbsParam.materialName()),
                0,
                coll5AbsorberOffsetInColl5,
                coll5Info,
                0,
                visible,
                G4Color::Gray(),
                solid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
    }
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
    bool visible             = config.getBool("muondegrader.visible",false);
    bool solid               = config.getBool("muondegrader.solid" ,false);
    int const verbosityLevel = config.getInt("muondegrader.verbosityLevel", 0);
    bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;

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
                      visible,
                      G4Color::Blue(),
                      solid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck
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
    bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;
    bool visible = config.getBool("pbar.visible",false);
    bool solid   = config.getBool("pbar.solid",false);

    // Throw exception if pbarwedge.build is used
    if ( config.hasName("pbarwedge.build") )
      {
        throw cet::exception("GEOM")<<
          " Variable pbarwedge.build is now deprecated. \n" <<
          " To use pbar wedge specify: pbar.Type = \"wedge\" \n" ;
      }


    PbarWindow const & pbarWindow = bl.getTS().getPbarWindow();
    G4Material* pbarMaterial  = findMaterialOrThrow( pbarWindow.material() );

    if ( pbarWindow.shape() == "disk" ) 
      {
        double pbarParams[5]  = { 0.0, pbarWindow.rOut(), pbarWindow.halfLength(), 0.0, CLHEP::twopi };
        
        VolumeInfo pbarInfo = nestTubs( "PbarAbs",
                                        pbarParams,
                                        pbarMaterial,
                                        0,
                                        pbarWindow.getLocal(),
                                        parent,
                                        0,
                                        visible,
                                        G4Color::Yellow(),
                                        solid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );
      }        
    else if( pbarWindow.shape() == "wedge" ) 
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
                      visible,
                      G4Color::Yellow(),
                      solid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck
                      );

      } 
    else if ( pbarWindow.shape() == "polycone" ) 
      {
        // Define polycone parameters
        vector<double> tmp_zPlanesDs3 {-0.50,-0.06,0.06,0.50};
        vector<double> tmp_rOuterDs3  (4,239.5);
        vector<double> tmp_rInnerDs3  {239.5,0.,0.,239.5};
        
        
        CLHEP::Hep3Vector polyPositionInMu2e( bl.getTS().getTS3_in().getGlobal() );
        
        nestPolycone( "PbarAbsPolycone",
                      PolyconsParams(tmp_zPlanesDs3,
                                     tmp_rInnerDs3,
                                     tmp_rOuterDs3 ),
                      pbarMaterial,
                      0,
                      polyPositionInMu2e - parent.centerInMu2e(),
                      parent,
                      0,
                      visible,
                      G4Colour::Yellow(),
                      solid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck
                      );
      }
    else 
      {
        throw cet::exception("GEOM")<<
          " Incorrect pbar window geometry requested! \n " ;
      }

  } // end Mu2eWorld::constructPbarWindow()

}

