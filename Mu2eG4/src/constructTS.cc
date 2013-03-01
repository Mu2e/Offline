//
// Free function to create Transport Solenoid
//
// $Id: constructTS.cc,v 1.12 2013/03/01 21:58:26 logash Exp $
// $Author: logash $
// $Date: 2013/03/01 21:58:26 $
//
// Original author KLG based on Mu2eWorld constructTS
//
// Notes:
// Construct the TS.  Parent volume is the air inside of the hall.

// C++ includes
#include <iostream>

// Framework includes
#include "cetlib/exception.h"

// Mu2e includes.
#include "Mu2eG4/inc/constructTS.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
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
//#include "G4CutTubs.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"

using namespace std;

namespace mu2e {

  void constructTS( VolumeInfo const & parent,
                    SimpleConfig const & _config
                    ){

    // Extract base parameters from config information.
    GeomHandle<Beamline> beamg;
    //double solenoidOffset = beamg->solenoidOffset();

    int const verbosityLevel = _config.getInt("coll.verbosityLevel", 0);

    double rTorus         = beamg->getTS().torusRadius();
    double rVac           = beamg->getTS().innerRadius();
    double rCryo          = beamg->getTS().outerRadius();
    double ts1HalfLength  = beamg->getTS().getTS1().getHalfLength();
    double ts3HalfLength  = beamg->getTS().getTS3().getHalfLength();
    double ts5HalfLength  = beamg->getTS().getTS5().getHalfLength();

    bool toyTSVisible = _config.getBool("toyTS.visible",true);
    bool toyTSSolid   = _config.getBool("toyTS.solid",true);

    double coll1HalfLength     = beamg->getTS().getColl1().getHalfLength();
    double coll31HalfLength    = beamg->getTS().getColl31().getHalfLength();
    double coll32HalfLength    = beamg->getTS().getColl32().getHalfLength();
    double coll5HalfLength     = beamg->getTS().getColl5().getHalfLength();

    double coll1InnerRadius1   = _config.getDouble("coll1.innerRadius1");
    double coll1InnerRadius2   = _config.getDouble("coll1.innerRadius2");
    double coll1InnerRadius3   = _config.getDouble("coll1.innerRadius3");

    double coll5InnerRadius    = _config.getDouble("coll5.innerRadius");
    double coll5MidRadius1     = _config.getDouble("coll5.midRadius1");
    double coll5MidRadius2     = _config.getDouble("coll5.midRadius2");
    double coll5OuterRadius    = _config.getDouble("coll5.outerRadius");

    double coll5HalfLengthU    = _config.getDouble("coll5.halfLengthU");
    double coll5HalfLengthD    = _config.getDouble("coll5.halfLengthD");

    MaterialFinder materialFinder(_config);

    G4Material* coll1Material1 = materialFinder.get("coll1.material1Name");
    G4Material* coll1Material2 = materialFinder.get("coll1.material2Name");
    G4Material* coll3Material  = materialFinder.get("coll3.materialName");
    // G4Material* coll5Material  = materialFinder.get("coll5.materialName");
    string coll5MaterialName         = _config.getString("coll5.materialName");
    string coll5AbsorberMaterialName = _config.getString("coll5.absorberMaterialName");

    // Special parameters for coll3
    double coll3RotationAngle    = _config.getDouble("coll3.rotationAngle");
    double coll3HoleRadius       = _config.getDouble("coll3.holeRadius");
    double coll3HoleHalfHeight   = _config.getDouble("coll3.holeHalfHeight");
    double coll3HoleDisplacement = _config.getDouble("coll3.holeDisplacement");

    bool collVisible         = _config.getBool("coll.visible",true);
    bool collSolid           = _config.getBool("coll.solid",true);
    bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;

    // For how all pieces are made from one of two types of material,
    // vacuum or average coils + cryostat material.
    G4Material* vacuumMaterial  = materialFinder.get("toyDS.insideMaterialName");
    G4Material* cryoMaterial    = materialFinder.get("toyDS.materialName");

    GeomHandle<VirtualDetector> vdg;
    double vdHalfLength = vdg->getHalfLength()*CLHEP::mm;

    // Computed quantities.
    //double ts5zOffset    = ( rTorus+ts5HalfLength);

    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

    // Build TS1.
    TubsParams ts1VacParams (   0.,  rVac, ts1HalfLength);
    TubsParams ts1CryoParams( rVac, rCryo, ts1HalfLength);

    VolumeInfo ts1VacInfo = nestTubs( "ToyTS1Vacuum",
                                      ts1VacParams,
                                      vacuumMaterial,
                                      beamg->getTS().getTS1().getRotation(),
                                      beamg->getTS().getTS1().getGlobal()-_hallOriginInMu2e,
                                      parent,
                                      0,
                                      toyTSVisible,
                                      G4Color::Red(),
                                      toyTSSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    VolumeInfo ts1CryoInfo = nestTubs( "ToyTS1Cryo",
                                       ts1CryoParams,
                                       cryoMaterial,
                                       beamg->getTS().getTS1().getRotation(),
                                       beamg->getTS().getTS1().getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       toyTSVisible,
                                       G4Color::Red(),
                                       toyTSSolid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    // Place collimator 1 (concentric cone which can be a cylinder when r1==r2
    // and a cylinder placed in ToyTS1Vacuum)

    // the cone (which can be a tube/cylinder) inside the outer tube/cylinder
    double coll1Param1[7] = { coll1InnerRadius1, coll1InnerRadius3,
                              coll1InnerRadius2, coll1InnerRadius3,
                              coll1HalfLength - 2.*vdHalfLength,
                              0.0, CLHEP::twopi };

    VolumeInfo coll1VacInfo = nestCons( "Coll11",
                                        coll1Param1,
                                        coll1Material1,
                                        0,
                                        beamg->getTS().getColl1().getLocal(),
                                        ts1VacInfo,
                                        0,
                                        collVisible,
                                        G4Color::Cyan(),
                                        collSolid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );

    TubsParams coll1Param2 ( coll1InnerRadius3,              rVac, coll1HalfLength-2.*vdHalfLength);
 
    VolumeInfo coll1Info2 = nestTubs( "Coll12",
                                      coll1Param2,
                                      coll1Material2,
                                      0,
                                      beamg->getTS().getColl1().getLocal(),
                                      ts1VacInfo,
                                      0,
                                      collVisible,
                                      G4Color::Blue(),
                                      collSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    // Build TS2.
    double ts2VacParams[5]  = { 0.0,   rVac, rTorus, 1.5*CLHEP::pi, CLHEP::halfpi };
    double ts2CryoParams[5] = { rVac, rCryo, rTorus, 1.5*CLHEP::pi, CLHEP::halfpi };

    // Position in the Mu2e coordintate system.
    G4ThreeVector ts2VacPosition( ts3HalfLength, 0., -rTorus);

    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
    G4RotationMatrix* ts2Rot = reg.add(G4RotationMatrix());
    ts2Rot->rotateX(90.0*CLHEP::degree);

    VolumeInfo ts2VacInfo = nestTorus("ToyTS2Vacuum",
                                      ts2VacParams,
                                      vacuumMaterial,
                                      ts2Rot,
                                      ts2VacPosition-_hallOriginInMu2e,
                                      parent,
                                      0,
                                      toyTSVisible,
                                      G4Color::Yellow(),
                                      toyTSSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    VolumeInfo ts2CryoInfo = nestTorus("ToyTS2Cryo",
                                       ts2CryoParams,
                                       cryoMaterial,
                                       ts2Rot,
                                       ts2VacPosition-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       toyTSVisible,
                                       G4Color::Yellow(),
                                       toyTSSolid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    // Build TS3.

    TubsParams ts3VacParams (   0.,  rVac, ts3HalfLength);
    TubsParams ts3CryoParams( rVac, rCryo, ts3HalfLength);

    VolumeInfo ts3VacInfo = nestTubs( "ToyTS3Vacuum",
                                      ts3VacParams,
                                      vacuumMaterial,
                                      beamg->getTS().getTS3().getRotation(),
                                      beamg->getTS().getTS3().getGlobal()-_hallOriginInMu2e,
                                      parent,
                                      0,
                                      toyTSVisible,
                                      G4Color::Green(),
                                      toyTSSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    VolumeInfo ts3CryoInfo = nestTubs( "ToyTS3Cryo",
                                       ts3CryoParams,
                                       cryoMaterial,
                                       beamg->getTS().getTS3().getRotation(),
                                       beamg->getTS().getTS3().getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       toyTSVisible,
                                       G4Color::Cyan(),
                                       toyTSSolid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    // Place collimator 3

    // Collimator 3 has peculiar shape, described in doc_db 853.
    // Construct this shape using boolean functions on solids

    // First, construct hole; make it slightly longer that any collimator
    double hDz = coll31HalfLength;
    if( hDz<coll32HalfLength ) hDz=coll32HalfLength;
    // Hole is the intersection of box and circles
    G4Box* coll3_hole_box = new G4Box("coll3_hole_box",
                                      coll3HoleRadius+5.0,coll3HoleHalfHeight,hDz+1.0);
    G4Tubs* coll3_hole_circle = new G4Tubs("coll3_hole_circle",
                                           0.0,coll3HoleRadius,hDz+1.0,
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
                                       0, rVac, coll31HalfLength-2.*vdHalfLength,
                                       0.0, CLHEP::twopi );

    G4Tubs* coll32_mother = new G4Tubs("Coll32_mother",
                                       0, rVac, coll32HalfLength-2.*vdHalfLength,
                                       0.0, CLHEP::twopi );

    coll31Info.solid = new G4SubtractionSolid(coll31Info.name,
                                              coll31_mother,
                                              coll3_hole,
                                              0,
                                              G4ThreeVector(0,coll3HoleDisplacement,0));

    coll32Info.solid = new G4SubtractionSolid(coll32Info.name,
                                              coll32_mother,
                                              coll3_hole,
                                              0,
                                              G4ThreeVector(0,coll3HoleDisplacement,0));

    // Now use finishNesting to place collimators 31 and 32

    G4RotationMatrix* coll31Rot = reg.add(G4RotationMatrix());
    G4RotationMatrix* coll32Rot = reg.add(G4RotationMatrix());
    coll31Rot->rotateZ(coll3RotationAngle*CLHEP::degree);
    coll32Rot->rotateZ(coll3RotationAngle*CLHEP::degree);

    finishNesting(coll31Info,
                  coll3Material,
                  coll31Rot,
                  beamg->getTS().getColl31().getLocal(),
                  ts3VacInfo.logical,
                  0,
                  collVisible,
                  G4Color::Gray(),
                  collSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);

    finishNesting(coll32Info,
                  coll3Material,
                  coll32Rot,
                  beamg->getTS().getColl32().getLocal(),
                  ts3VacInfo.logical,
                  0,
                  collVisible,
                  G4Color::Gray(),
                  collSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck);

    // Place Pbar absorber between Coll31 and Coll32
    // Pbar absorber is made of two pieces:
    //  -- vacuum wall, which covers the whole inner part of TS3
    //     it is controlled by pbar.* parameters
    //  -- wedge, which starts near center and extends upward
    //     it is controlled by pbarwedge.* parameters

    // -- vacuum wall
    
    double pbarHalfLength     = _config.getDouble("pbar.halfLength");
    G4Material* pbarMaterial  = materialFinder.get("pbar.materialName");
    double pbarParams[5]  = { 0.0,   rVac, pbarHalfLength, 0.0, CLHEP::twopi };

    VolumeInfo pbarInfo = nestTubs( "PbarAbs",
                                    pbarParams,
                                    pbarMaterial,
                                    0,
                                    G4ThreeVector(0.,0.,0.),
                                    ts3VacInfo,
                                    0,
                                    collVisible,
                                    G4Color::Yellow(),
                                    collSolid,
                                    forceAuxEdgeVisible,
                                    placePV,
                                    doSurfaceCheck
                                    );

    // -- pbar wedge
    
    bool addPbarWedge    = _config.getBool("pbarwedge.build",false);
    double pbarWedge_y0  = _config.getDouble("pbarwedge.y0",0);
    double pbarWedge_y1  = _config.getDouble("pbarwedge.y1",rVac-1);
    double pbarWedge_dz0 = _config.getDouble("pbarwedge.dz0",0.001);
    double pbarWedge_dz1 = _config.getDouble("pbarwedge.dz1",0.001);

    if( addPbarWedge ) {

      VolumeInfo pbarWedgeInfo;
      
      pbarWedgeInfo.name = "PbarAbsWedge";

      double pbarWedge_dz = ( pbarWedge_dz0<pbarWedge_dz1 ) ? pbarWedge_dz1 : pbarWedge_dz0;
      double pbarWedge_h = pbarWedge_y1 - pbarWedge_y0;
      
      double pbarWedge_dy = (pbarWedge_y1 + pbarWedge_y0)/2.;
      
      G4Tubs *pbarWedge_disk = new G4Tubs("PbarAbsWedge_disk",
					  0,rVac,pbarWedge_dz/2.,0,CLHEP::twopi);

      G4Trd *pbarWedge_trd = new G4Trd("PbarAbsWedge_trd",
				       rVac,rVac,
				       pbarWedge_dz0/2.,pbarWedge_dz1/2.,
				       pbarWedge_h/2.);
      
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
		    G4ThreeVector(0.,0.,pbarWedge_dz/2+pbarHalfLength),
		    ts3VacInfo.logical,
		    0,
		    collVisible,
		    G4Color::Yellow(),
		    collSolid,
		    forceAuxEdgeVisible,
		    placePV,
		    doSurfaceCheck
		    );

    }

    // Build TS4.

    double ts4VacParams[5]  = { 0.0,  rVac, rTorus, CLHEP::halfpi, CLHEP::halfpi };
    double ts4CryoParams[5] = {rVac, rCryo, rTorus, CLHEP::halfpi, CLHEP::halfpi };

    // Position in the Mu2e coordintate system.

    G4ThreeVector ts4VacPosition( -ts3HalfLength, 0., rTorus);

    G4RotationMatrix* ts4Rot = reg.add(G4RotationMatrix());
    ts4Rot->rotateX(90.0*CLHEP::degree);

    VolumeInfo ts4VacInfo = nestTorus("ToyTS4Vacuum",
                                      ts4VacParams,
                                      vacuumMaterial,
                                      ts4Rot,
                                      ts4VacPosition-_hallOriginInMu2e,
                                      parent,
                                      0,
                                      toyTSVisible,
                                      G4Color::Yellow(),
                                      toyTSSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    VolumeInfo ts4CryoInfo = nestTorus("ToyTS4Cryo",
                                       ts4CryoParams,
                                       cryoMaterial,
                                       ts4Rot,
                                       ts4VacPosition-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       toyTSVisible,
                                       G4Color::Yellow(),
                                       toyTSSolid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    // Build TS5.

    TubsParams ts5VacParams (   0.,  rVac, ts5HalfLength);
    TubsParams ts5CryoParams( rVac, rCryo, ts5HalfLength);

    VolumeInfo ts5VacInfo = nestTubs( "ToyTS5Vacuum",
                                      ts5VacParams,
                                      vacuumMaterial,
                                      beamg->getTS().getTS5().getRotation(),
                                      beamg->getTS().getTS5().getGlobal()-_hallOriginInMu2e,
                                      parent,
                                      0,
                                      toyTSVisible,
                                      G4Color::Red(),
                                      toyTSSolid,
                                      forceAuxEdgeVisible,
                                      placePV,
                                      doSurfaceCheck
                                      );

    VolumeInfo ts5CryoInfo = nestTubs( "ToyTS5Cryo",
                                       ts5CryoParams,
                                       cryoMaterial,
                                       beamg->getTS().getTS5().getRotation(),
                                       beamg->getTS().getTS5().getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       toyTSVisible,
                                       G4Color::Red(),
                                       toyTSSolid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    // Place collimator 5

    if ( coll5OuterRadius > rVac || 
         coll5MidRadius2  > coll5OuterRadius - 2.*vdHalfLength || 
         coll5MidRadius1  >= coll5MidRadius2 || 
         coll5InnerRadius > coll5MidRadius1 ||
         coll5InnerRadius < 0.0 ) {

      throw cet::exception("GEOM")<< " constructTS: wrong coll5 radii: " 
                                  << "\n rVac             : " << rVac
                                  << "\n coll5OuterRadius : " << coll5OuterRadius
                                  << "\n coll5MidRadius2  : " << coll5MidRadius2
                                  << "\n coll5MidRadius1  : " << coll5MidRadius1
                                  << "\n coll5InnerRadiu  : " << coll5InnerRadius
                                  << "\n";
    }

    if ( coll5HalfLengthU < 0.0 || coll5HalfLengthD < 0.0 ||
         coll5HalfLengthU + coll5HalfLengthD > coll5HalfLength - 2.*vdHalfLength) {

      throw cet::exception("GEOM")<< " constructTS: wrong coll5 longitudinal params " 
                                  << "\n coll5HalfLength   : " << coll5HalfLength
                                  << "\n coll5HalfLengthU  : " << coll5HalfLengthU
                                  << "\n coll5HalfLengthD  : " << coll5HalfLengthD
                                  << "\n";
    }

    if ( verbosityLevel > 0) {

      cout << __func__ << " TS5  OffsetInMu2e    : "
           << beamg->getTS().getTS5().getGlobal() << endl;
      cout << __func__ << " Coll5 local offset   : " 
           << beamg->getTS().getColl5().getLocal() << endl;
      cout << __func__ << " beamg->getTS().getTS5().getRotation(): "
           << beamg->getTS().getTS5().getRotation() << endl;

    }

    CLHEP::Hep3Vector coll5OffsetInMu2e = beamg->getTS().getTS5().getGlobal() + 
      ( ( beamg->getTS().getTS5().getRotation() != 0x0 ) ?
        *(beamg->getTS().getTS5().getRotation()) * beamg->getTS().getColl5().getLocal() : 
        beamg->getTS().getColl5().getLocal() );

    if ( verbosityLevel > 0) {

      cout << __func__ << "  coll5OffsetInMu2e    : "
           << coll5OffsetInMu2e << endl;
      cout << __func__ << "  Coll5 calc local offset : "
           << coll5OffsetInMu2e - beamg->getTS().getTS5().getGlobal() << endl;

    }

    // the most outer part (with Virtual Detectors on the outer surfaces of the Coll5)

    Tube coll5Param(coll5MaterialName,
                    coll5OffsetInMu2e,
                    coll5InnerRadius,
                    coll5OuterRadius - 2.*vdHalfLength,
                    coll5HalfLength - 2.*vdHalfLength);

    VolumeInfo coll5Info = nestTubs( "Coll5",
                                     coll5Param.getTubsParams(),
                                     findMaterialOrThrow(coll5Param.materialName()),
                                     0,
                                     coll5Param.originInMu2e() - ts5VacInfo.centerInMu2e(),
                                     ts5VacInfo,
                                     0,
                                     collVisible,
                                     G4Color::Blue(),
                                     collSolid,
                                     forceAuxEdgeVisible,
                                     placePV,
                                     doSurfaceCheck
                                     );

    // decide if we need to add absorber parts

    if ( coll5HalfLengthD > 0.0 ) {
      
      // + downstream absorber
      
      CLHEP::Hep3Vector coll5AbsorberOffsetInColl5 = 
        CLHEP::Hep3Vector (0.0, 0.0, coll5HalfLength - 2.*vdHalfLength - coll5HalfLengthD);

      CLHEP::Hep3Vector coll5AbsorberOffsetInMu2e = coll5OffsetInMu2e +
        ( ( beamg->getTS().getTS5().getRotation() != 0x0 ) ?
          *(beamg->getTS().getTS5().getRotation()) * coll5AbsorberOffsetInColl5 :
          coll5AbsorberOffsetInColl5 );

      Tube coll5DAbsParam(coll5AbsorberMaterialName,
                          coll5AbsorberOffsetInMu2e,
                          coll5InnerRadius,
                          coll5OuterRadius - 2.*vdHalfLength,
                          coll5HalfLengthD);

      VolumeInfo coll5DAbsorberInfo = nestTubs( "coll5DAbsorber",
                                                coll5DAbsParam.getTubsParams(),
                                                findMaterialOrThrow(coll5DAbsParam.materialName()),
                                                0,
                                                coll5AbsorberOffsetInColl5,
                                                coll5Info,
                                                0,
                                                collVisible,
                                                G4Color::Gray(),
                                                collSolid,
                                                forceAuxEdgeVisible,
                                                placePV,
                                                doSurfaceCheck
                                                );
    }

    if ( coll5HalfLengthU > 0.0 ) {

      // + upstream absorber
    
      CLHEP::Hep3Vector coll5AbsorberOffsetInColl5 = 
        CLHEP::Hep3Vector (0.0, 0.0, - coll5HalfLength + 2.*vdHalfLength + coll5HalfLengthU);

      CLHEP::Hep3Vector coll5AbsorberOffsetInMu2e = coll5OffsetInMu2e +
        ( ( beamg->getTS().getTS5().getRotation() != 0x0 ) ?
          *(beamg->getTS().getTS5().getRotation()) * coll5AbsorberOffsetInColl5 :
          coll5AbsorberOffsetInColl5 );

      Tube coll5AbsParam(coll5AbsorberMaterialName,
                          coll5AbsorberOffsetInMu2e,
                          coll5InnerRadius,
                          coll5OuterRadius - 2.*vdHalfLength,
                          coll5HalfLengthU);

      VolumeInfo coll5UAbsorberInfo = nestTubs( "coll5UAbsorber",
                                                coll5AbsParam.getTubsParams(),
                                                findMaterialOrThrow(coll5AbsParam.materialName()),
                                                0,
                                                coll5AbsorberOffsetInColl5,
                                                coll5Info,
                                                0,
                                                collVisible,
                                                G4Color::Gray(),
                                                collSolid,
                                                forceAuxEdgeVisible,
                                                placePV,
                                                doSurfaceCheck
                                                );
    }

    if ( coll5MidRadius1 > coll5InnerRadius ) {

      // + "inner radius" absorber
    
      CLHEP::Hep3Vector coll5AbsorberOffsetInColl5 =
        CLHEP::Hep3Vector (0.0, 0.0,
                           coll5HalfLengthU - coll5HalfLengthD);

      CLHEP::Hep3Vector coll5AbsorberOffsetInMu2e = coll5OffsetInMu2e +
        ( ( beamg->getTS().getTS5().getRotation() != 0x0 ) ?
          *(beamg->getTS().getTS5().getRotation()) * coll5AbsorberOffsetInColl5 :
          coll5AbsorberOffsetInColl5 );

      Tube coll5AbsParam(coll5AbsorberMaterialName,
                         coll5AbsorberOffsetInMu2e,
                         coll5InnerRadius,
                         coll5MidRadius1,
                         coll5HalfLength - 2.*vdHalfLength - 
                         coll5HalfLengthD - coll5HalfLengthU);

      VolumeInfo coll5IAbsorberInfo = nestTubs( "coll5IAbsorber",
                                                coll5AbsParam.getTubsParams(),
                                                findMaterialOrThrow(coll5AbsParam.materialName()),
                                                0,
                                                coll5AbsorberOffsetInColl5,
                                                coll5Info,
                                                0,
                                                collVisible,
                                                G4Color::Gray(),
                                                collSolid,
                                                forceAuxEdgeVisible,
                                                placePV,
                                                doSurfaceCheck
                                                );
    }

    if ( coll5MidRadius2 < coll5OuterRadius - 2.*vdHalfLength) {

      // + "outer radius" absorber
    
      CLHEP::Hep3Vector coll5AbsorberOffsetInColl5 =
        CLHEP::Hep3Vector (0.0, 0.0,
                           coll5HalfLengthU - coll5HalfLengthD);

      CLHEP::Hep3Vector coll5AbsorberOffsetInMu2e = coll5OffsetInMu2e +
        ( ( beamg->getTS().getTS5().getRotation() != 0x0 ) ?
          *(beamg->getTS().getTS5().getRotation()) * coll5AbsorberOffsetInColl5 :
          coll5AbsorberOffsetInColl5 );

      Tube coll5AbsParam(coll5AbsorberMaterialName,
                         coll5AbsorberOffsetInMu2e,
                         coll5MidRadius2,
                         coll5OuterRadius - 2.*vdHalfLength,
                         coll5HalfLength - 2.*vdHalfLength - 
                         coll5HalfLengthD - coll5HalfLengthU);

      VolumeInfo coll5OAbsorberInfo = nestTubs( "coll5OAbsorber",
                                                coll5AbsParam.getTubsParams(),
                                                findMaterialOrThrow(coll5AbsParam.materialName()),
                                                0,
                                                coll5AbsorberOffsetInColl5,
                                                coll5Info,
                                                0,
                                                collVisible,
                                                G4Color::Gray(),
                                                collSolid,
                                                forceAuxEdgeVisible,
                                                placePV,
                                                doSurfaceCheck
                                                );
    }

    // Add muon degrader to the center of TS5
    // Degrader is made of several layers. 
    // Each layer is intersection of cylinder and trapezoid.
    
    bool addDegrader  = _config.getBool("muondegrader.build",false);
    vector<double> degraderR, degraderDZB, degraderDZT, degraderPhi;
    _config.getVectorDouble("muondegrader.R", degraderR, vector<double>() );
    _config.getVectorDouble("muondegrader.DZB", degraderDZB, vector<double>() );
    _config.getVectorDouble("muondegrader.DZT", degraderDZT, vector<double>() );
    _config.getVectorDouble("muondegrader.Phi", degraderPhi, vector<double>() );

    if( degraderR.size()!=degraderDZB.size() || degraderR.size()!=degraderDZT.size() ||
	degraderR.size()!=degraderPhi.size() ) {
      cout << __func__ << " Warning: MuonDegrader is not build - dimensions don't match." << endl;
      addDegrader = false;
    }

    if( addDegrader ) {

      G4Material* degraderMaterial  = materialFinder.get("muondegrader.materialName");
      
      for( unsigned int i=0; i<degraderR.size(); ++i ) {
	
	VolumeInfo degraderInfo;
	
	ostringstream dname;  dname << "MuonDegrader" << i;
	degraderInfo.name = dname.str();

	ostringstream dsname1;  dsname1 << "MuonDegrader_disk" << i;
	ostringstream dsname2;  dsname2 << "MuonDegrader_trd" << i;

	double R1 = degraderR[i];
	double R2 = ( i<=0 ) ? coll5InnerRadius : degraderR[i-1];
	G4Tubs *degrader_disk = new G4Tubs(dsname1.str(),R1,R2,
					   degraderDZB[i]/2.0,
					   (1.5-degraderPhi[i])*CLHEP::pi, 
					   degraderPhi[i]*CLHEP::twopi);

	G4Trd *degrader_trd = new G4Trd(dsname2.str(),
					coll5InnerRadius,coll5InnerRadius,
					degraderDZT[i]/2.0,degraderDZB[i]/2.0,
					R2/2);
	
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
		      coll5Param.originInMu2e() - ts5VacInfo.centerInMu2e(),
		      ts5VacInfo.logical,
		      0,
		      collVisible,
		      G4Color::Blue(),
		      collSolid,
		      forceAuxEdgeVisible,
		      placePV,
		      doSurfaceCheck
		      );
      }

    }
    
  } // end Mu2eWorld::constructTS

}
