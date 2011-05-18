//
// Free function to create Transport Solenoid
//
// $Id: constructTS.cc,v 1.4 2011/05/18 14:21:44 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/18 14:21:44 $
//
// Original author KLG based on Mu2eWorld constructTS
//
// Notes:
// Construct the TS.  Parent volume is the air inside of the hall.

// C++ includes
//#include <iostream>

// Mu2e includes.
#include "Mu2eG4/inc/constructTS.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestTorus.hh"
#include "Mu2eG4/inc/nestCons.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "VirtualDetectorGeom/inc/VirtualDetector.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"


using namespace std;

namespace mu2e {

  void constructTS( const VolumeInfo& parent,
                    SimpleConfig const * const _config
                    ){

    // Extract base parameters from config information.
    GeomHandle<Beamline> beamg;
    //double solenoidOffset = beamg->solenoidOffset();
    double rTorus         = beamg->getTS().torusRadius();
    double rVac           = beamg->getTS().innerRadius();
    double rCryo          = beamg->getTS().outerRadius();
    double ts1HalfLength  = beamg->getTS().getTS1().getHalfLength();
    double ts3HalfLength  = beamg->getTS().getTS3().getHalfLength();
    double ts5HalfLength  = beamg->getTS().getTS5().getHalfLength();

    bool toyTSVisible = _config->getBool("toyTS.visible",true);
    bool toyTSSolid   = _config->getBool("toyTS.solid",true);

    double coll1HalfLength     = beamg->getTS().getColl1().getHalfLength();
    double coll31HalfLength    = beamg->getTS().getColl31().getHalfLength();
    double coll32HalfLength    = beamg->getTS().getColl32().getHalfLength();
    double coll5HalfLength     = beamg->getTS().getColl5().getHalfLength();

    double coll1InnerRadius1   = _config->getDouble("coll1.innerRadius1");
    double coll1InnerRadius2   = _config->getDouble("coll1.innerRadius2");
    double coll5InnerRadius    = _config->getDouble("coll5.innerRadius");

    MaterialFinder materialFinder(*_config);

    G4Material* coll1Material  = materialFinder.get("coll1.materialName");
    G4Material* coll3Material  = materialFinder.get("coll3.materialName");
    G4Material* coll5Material  = materialFinder.get("coll5.materialName");

    // Special parameters for coll3
    double coll3RotationAngle    = _config->getDouble("coll3.rotationAngle");
    double coll3HoleRadius       = _config->getDouble("coll3.holeRadius");
    double coll3HoleHalfHeight   = _config->getDouble("coll3.holeHalfHeight");
    double coll3HoleDisplacement = _config->getDouble("coll3.holeDisplacement");

    bool collVisible         = _config->getBool("coll.visible",true);
    bool collSolid           = _config->getBool("coll.solid",true);
    bool forceAuxEdgeVisible = _config->getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config->getBool("g4.doSurfaceCheck",false);
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

    // Place collimator 1

    double coll1Param[7] = { coll1InnerRadius1, rVac,
			     coll1InnerRadius2, rVac,
			     coll1HalfLength-2*vdHalfLength,
			     0.0, 360.0*CLHEP::degree };

    VolumeInfo coll1VacInfo = nestCons( "Coll1",
                                        coll1Param,
                                        coll1Material,
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
    double ts2VacParams[5]  = { 0.0,   rVac, rTorus, 1.5*M_PI, 0.5*M_PI };
    double ts2CryoParams[5] = { rVac, rCryo, rTorus, 1.5*M_PI, 0.5*M_PI };

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
					   0.0, 360.0*CLHEP::degree );
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
				       0, rVac, coll31HalfLength-2*vdHalfLength,
				       0.0, 360.0*CLHEP::degree );

    G4Tubs* coll32_mother = new G4Tubs("Coll32_mother",
				       0, rVac, coll32HalfLength-2*vdHalfLength,
				       0.0, 360.0*CLHEP::degree );

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

    double pbarHalfLength     = _config->getDouble("pbar.halfLength");
    G4Material* pbarMaterial  = materialFinder.get("pbar.materialName");
    double pbarParams[5]  = { 0.0,   rVac, pbarHalfLength, 0.0, 360.0*CLHEP::degree };

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

    // Build TS4.

    double ts4VacParams[5]  = { 0.0,  rVac, rTorus, 0.5*M_PI, 0.5*M_PI };
    double ts4CryoParams[5] = {rVac, rCryo, rTorus, 0.5*M_PI, 0.5*M_PI };

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

    double coll5Param[5] = { coll5InnerRadius, rVac,
			     coll5HalfLength-2*vdHalfLength, 0.0, 360.0*CLHEP::degree };

    VolumeInfo coll5VacInfo = nestTubs( "Coll5",
                                        coll5Param,
                                        coll5Material,
                                        0,
                                        beamg->getTS().getColl5().getLocal(),
                                        ts5VacInfo,
                                        0,
                                        collVisible,
                                        G4Color::Blue(),
                                        collSolid,
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );

  } // end Mu2eWorld::constructTS

}
