//
// Construct the Mu2e G4 world and serve information about that world.
//
// $Id: Mu2eWorld.cc,v 1.56 2010/09/29 19:37:58 logash Exp $
// $Author: logash $ 
// $Date: 2010/09/29 19:37:58 $
//
// Original author Rob Kutschke
//
//  Heirarchy is:
//  0      World (air)
//  1      Earthen Overburden
//  2      Concrete walls of the hall
//  3      Air inside the hall
//  4      Effective volume representing the DS coils+cryostats.
//  4      DS Vaccum
//
//  4      Effective volume representing the PS coils+cryostats.
//  4      PS Vacuum
//
//  The Earth overburden is modeled in two parts: a box that extends
//  to the surface of the earth plus a cap above grade.  The cap is shaped
//  as a G4Paraboloid.
//
// Notes:
// 1) When a G4VisAttributes is given to a G4LogicalVolume, the logical volume
//    object does not take ownership of the G4VisAttributes object; so the caller
//    who called SetVisAttributes needs to manage the lifetime of the G4VisAttributes
//    object.  We do this by putting them into a list that has a lifetime as long
//    as Mu2eWorld.  Note that push_back() on a list does not invalidate pointers
//    to previous entries in the list.
//

// C++ includes
#include <iostream>
#include <vector>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/StrawPlacer.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/VirtualDetectorSD.hh"
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/ITrackerBuilder.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "TargetGeom/inc/Target.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/TransportSolenoid.hh"
#include "VirtualDetectorGeom/inc/VirtualDetector.hh"
#include "Mu2eG4/inc/constructLTracker.hh"
#include "Mu2eG4/inc/constructTTracker.hh"
#include "Mu2eG4/inc/constructDummyTracker.hh"
#include "Mu2eG4/inc/constructStoppingTarget.hh"
#include "Mu2eG4/inc/constructDummyStoppingTarget.hh"
#include "Mu2eG4/inc/constructCalorimeter.hh"

// G4 includes
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Paraboloid.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Torus.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"
#include "G4ClassicalRK4.hh"
#include "G4ImplicitEuler.hh"
#include "G4ExplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"

#include "Mu2eG4/inc/DSField.hh"
#include "Mu2eG4/inc/FieldMgr.hh"

using namespace std;

namespace mu2e {

  Mu2eWorld::Mu2eWorld():
    _cosmicReferencePoint(),
    _mu2eOrigin(),
    _info(){
  }
  
  Mu2eWorld::~Mu2eWorld(){

    // Do not destruct the solids, logical volumes or physical volumes.
    // G4 looks after that itself.

  }

  // A helper function for debugging.  Print a subset of the physical volume store.
  void printPhys() {
    G4PhysicalVolumeStore* pstore = G4PhysicalVolumeStore::GetInstance();
    int n(0);
    for ( std::vector<G4VPhysicalVolume*>::const_iterator i=pstore->begin(); i!=pstore->end(); i++){
      cout << "Physical Volume: "
           << setw(5) << n++
           << (*i)->GetName()
           << endl;
      if ( n > 25 ) break;
    }

  }


  // This is the callback called by G4.
  WorldInfo const* Mu2eWorld::construct(){

    // In case this is called a second time within a job.
    _volumeInfoList.clear();
    _visAttributes.clear();
    
    // Get access to the master geometry system and its run time config.
    edm::Service<GeometryService> geom;
    SimpleConfig const& config = geom->config();
    _config = &config;

    // Construct all of the Mu2e world, hall, detectors, beamline ...
    constructWorld();

    return &_info;
  }

  // Construct all of the Mu2e world, hall, detectors, beamline ...
  void Mu2eWorld::constructWorld(){

    // If you play with the order of these calls, you may break things.
    defineMu2eOrigin();
    VolumeInfo::setMu2eOriginInWorld( _mu2eOrigin );
    VirtualDetectorSD::setMu2eOriginInWorld( _mu2eOrigin );
    CaloCrystalSD::setMu2eOriginInWorld( _mu2eOrigin );
    VolumeInfo dirtInfo  = constructDirt();
    VolumeInfo hallInfo  = constructHall( dirtInfo );
    constructDS(hallInfo);
    constructTS(hallInfo);
    constructPS(hallInfo);
    constructVD();
    VolumeInfo trackerInfo = constructTracker();
    VolumeInfo targetInfo  = constructTarget();
    constructProtonAbs();  //this is to construct proton absorber

    // These are just placeholders for now - and might be misnamed.
    constructCal();
    constructMagnetYoke();
    constructCRV();

    // This does real work.
    constructSteel(hallInfo);

    // Hack alert: These belong in constructTracker and constructTarget.
    trackerInfo.name = "TrackerMother";
    targetInfo.name  = "StoppingTargetMother";
    addVolInfo(trackerInfo);
    addVolInfo(targetInfo);

    edm::LogInfo log("GEOM");
    log << "Mu2e Origin:          " << _mu2eOrigin           << "\n";
    log << "Mu2e Detector Origin: " << _mu2eDetectorOrigin   << "\n";
    log << "Cosmic Ref:           " << _cosmicReferencePoint << "\n";
    log << "Hall Origin in Mu2e:  " << _hallOriginInMu2e     << "\n";

    // Create magnetic fields and managers only after all volumes have been defined.
    // Eventually constructBFieldAndManagers2() will be the only BField constructor.
    if( _config->getString("bfield.format","GMC") == "GMC" ) {
      constructBFieldAndManagers();
    } else {
      constructBFieldAndManagers2();
    }
    constructStepLimiters();

  }

  // Convert to base units for all of the items in the vector.
  void Mu2eWorld::setUnits( vector<double>& V, G4double unit ){
    for ( vector<double>::iterator b=V.begin(), e=V.end();
          b!=e; ++b){
      *b *= unit;
    }
  }

  void Mu2eWorld::defineMu2eOrigin(){

    // Dimensions of the world.
    vector<double> worldHLen;
    _config->getVectorDouble("world.halfLengths", worldHLen, 3);

    // Floor thickness.
    double floorThick = _config->getDouble("hall.floorThick");

    // Top of the floor in G4 world coordinates.
    double yFloor = -worldHLen[1] + floorThick;

    // The height above the floor of the y origin of the Mu2e coordinate system.
    double yOriginHeight = _config->getDouble("world.mu2eOrigin.height" )*CLHEP::mm;

    // Position of the origin of the mu2e coordinate system
    _mu2eOrigin = G4ThreeVector( 
                                _config->getDouble("world.mu2eOrigin.xoffset")*CLHEP::mm,
                                yFloor + yOriginHeight,
                                _config->getDouble("world.mu2eOrigin.zoffset")*CLHEP::mm
                                );

    // Origin used to construct the MECO detector.
    // Magic number to fix:
    _mu2eDetectorOrigin = _mu2eOrigin + G4ThreeVector( -3904., 0., 12000.);

    double ceilingThick     = _config->getDouble("hall.ceilingThick");
    double overburdenDepth  = _config->getDouble("dirt.overburdenDepth");
    double capHalfHeight    = _config->getDouble("dirt.capHalfHeight");
    vector<double> hallInHLen;
    _config->getVectorDouble("hall.insideHalfLengths",hallInHLen,3);

    // Bottom of the ceiling in G4 world coordinates.
    double yCeilingInSide = yFloor + 2.*hallInHLen[1];
    
    // Top of the ceiling in G4 world coordinates.
    double yCeilingOutside  = yCeilingInSide + ceilingThick;

    // Surface of the earth in G4 world coordinates.
    double ySurface  = yCeilingOutside + overburdenDepth;

    // Top of the world.
    double yEverest = ySurface + 2.*capHalfHeight;

    // Build the reference points that others will use.
    _cosmicReferencePoint = G4ThreeVector( 0., yEverest, 0.);

    // Selfconsistency check.
    if ( yEverest > 2.*worldHLen[1] ){
      throw cms::Exception("GEOM")
        << "Top of the world is outside of the world volume! \n";
    }

  }  // end of Mu2eWorld::defineMu2eOrigin


  // Construct hall walls and hall interior inside the earthen overburden.
  VolumeInfo Mu2eWorld::constructHall ( const VolumeInfo& parent ){

    // A helper class.
    MaterialFinder materialFinder(*_config);

    // Dimensions of the world.
    vector<double> worldHLen;
    _config->getVectorDouble("world.halfLengths", worldHLen, 3);

    vector<double> hallInHLen;
    _config->getVectorDouble("hall.insideHalfLengths",hallInHLen,3);

    // Floor thickness.
    double ceilingThick = _config->getDouble("hall.ceilingThick");
    double floorThick   = _config->getDouble("hall.floorThick");
    double wallThick    = _config->getDouble("hall.wallThick");

    // Top of the floor in G4 world coordinates.
    double yFloor = -worldHLen[1] + floorThick;

    // Position of the center of the hall in the world volume.
    vector<double> hallPosition;
    _config->getVectorDouble("hall.offset", hallPosition,3);
    double hallY0 = yFloor + hallInHLen[1] + hallPosition[1];

    // Materials for the hall walls and the interior of the hall
    G4Material* wallMaterial = materialFinder.get("hall.wallMaterialName");
    G4Material* hallMaterial = materialFinder.get("hall.insideMaterialName");

    // Half lengths of the exterior of the concrete for the hall walls.
    double hallOutHLen[3] ={
      hallInHLen[0] + wallThick,
      hallInHLen[1] + ( ceilingThick + floorThick )/2.,
      hallInHLen[2] + wallThick
    };
    
    // Center of the concrete volume in the coordinate system of the dirt.
    G4ThreeVector wallOffset = 
      G4ThreeVector(hallPosition[0], hallY0, hallPosition[2]) - parent.centerInParent;

    // Origin of the hall air volume in the system of the hall concrete volume.
    G4ThreeVector hallOffset( 0., (floorThick-ceilingThick)/2., 0.);

    bool hallVisible = _config->getBool("hall.visible",true);
    bool hallSolid   = _config->getBool("hall.solid",false);

    // Concrete walls of the hall.
    VolumeInfo wallInfo = nestBox( "HallWalls",
                                   hallOutHLen,
                                   wallMaterial,
                                   0,
                                   wallOffset,
                                   parent,
                                   0,
                                   hallVisible,
                                   G4Colour::Red(),
                                   hallSolid
                                   );
    
    // Air volume inside of the hall.
    VolumeInfo hallInfo = nestBox( "HallAir",
                                   hallInHLen,
                                   hallMaterial,
                                   0,
                                   hallOffset,
                                   wallInfo,
                                   0,
                                   hallVisible,
                                   G4Colour::Red(),
                                   hallSolid
                                   );

    // Define the hall origin in Mu2e coordinates.
    _hallOriginInMu2e = hallInfo.centerInWorld - _mu2eOrigin;

    return hallInfo;
  }  // end of Mu2eWorld::constructHall

  // Construct world mother volume and the earthen overburden.
  VolumeInfo Mu2eWorld::constructDirt (){

    // A helper class.
    MaterialFinder materialFinder(*_config);
    
    // Dimensions and material of the world.
    vector<double> worldHLen;
    _config->getVectorDouble("world.halfLengths", worldHLen, 3);
    G4Material* worldMaterial = materialFinder.get("world.materialName");

    const string worldName("World");

    bool worldBoxVisible = _config->getBool("world.boxVisible",true);
    bool worldBoxSolid   = _config->getBool("world.boxSolid",false);

    // Construct the world volume.  The dummy is needed because the interface
    // to nestBox requires a mother even if this is the top level.
    VolumeInfo dummy("MotherofTheWorld", G4ThreeVector(), G4ThreeVector());
    VolumeInfo worldInfo = nestBox (worldName, 
                                    worldHLen, 
                                    worldMaterial, 
                                    0, 
                                    G4ThreeVector(), 
                                    dummy,
                                    0,
                                    worldBoxVisible,
                                    worldBoxSolid
                                    );
    _info.worldPhys  = worldInfo.physical;

    // Get parameters related to the overall dimensions of the hall and to
    // the earthen overburden.
    double floorThick           = CLHEP::mm * _config->getDouble("hall.floorThick");
    double ceilingThick         = CLHEP::mm * _config->getDouble("hall.ceilingThick");
    double wallThick            = CLHEP::mm * _config->getDouble("hall.wallThick");
    double overburdenDepth      = CLHEP::mm * _config->getDouble("dirt.overburdenDepth");
    vector<double> hallInHLen;
    _config->getVectorDouble("hall.insideHalfLengths",hallInHLen,3);

    // Derived parameters.
    G4Material* dirtMaterial = materialFinder.get("dirt.overburdenMaterialName");

    // Top of the floor in G4 world coordinates.
    double yFloor   = -worldHLen[1] + floorThick;

    // The height above the floor of the y origin of the Mu2e coordinate system.
    //    double yOriginHeight = _config->getDouble("world.mu2eOrigin.height" )*CLHEP::mm;

    // Bottom of the ceiling in G4 world coordinates.
    double yCeilingInSide = yFloor + 2.*hallInHLen[1];
    
    // Top of the ceiling in G4 world coordinates.
    double yCeilingOutside  = yCeilingInSide + ceilingThick;

    // Surface of the earth in G4 world coordinates.
    double ySurface  = yCeilingOutside + overburdenDepth;
    
    // Half length and y origin of the dirt box.
    double yLDirt = ( ySurface + worldHLen[1] )/2.;
    double y0Dirt = -worldHLen[1] + yLDirt;
    
    // Center of the dirt box, in the G4 world system.
    G4ThreeVector dirtOffset(0.,y0Dirt,0.);
    
    // Half lengths of the dirt box.
    double dirtHLen[3] = { worldHLen[0], yLDirt, worldHLen[2] };

    bool dirtVisible    = _config->getBool("dirt.visible",true);
    bool dirtSolid      = _config->getBool("dirt.solid",false);
    bool dirtCapVisible = _config->getBool("dirt.capVisible",true);
    bool dirtCapSolid   = _config->getBool("dirt.capSolid",false);

    // Main body of dirt around the hall.
    VolumeInfo dirtInfo = nestBox( "DirtBody",
                                   dirtHLen,
                                   dirtMaterial,
                                   0,
                                   dirtOffset,
                                   worldInfo,
                                   0,
                                   dirtVisible,
                                   G4Colour::Magenta(),
                                   dirtSolid
                                   );

    // Dirt cap is modeled as a paraboloid.
    double capHalfHeight = _config->getDouble("dirt.capHalfHeight");
    double capBottomR    = _config->getDouble("dirt.capBottomRadius");
    double capTopR       = _config->getDouble("dirt.capTopRadius");

    double dsz0          = _config->getDouble("toyDS.z0");

    GeomHandle<Beamline> beamg;
    double solenoidOffset = beamg->solenoidOffset();

    G4ThreeVector dirtCapOffset( -solenoidOffset, ySurface+capHalfHeight, dsz0+_mu2eOrigin.z());

    G4RotationMatrix* dirtCapRot = new G4RotationMatrix();
    dirtCapRot->rotateX( -90*CLHEP::degree);

    string dirtCapName("DirtCap");

    // Construct the cap.
    VolumeInfo dirtCapInfo( dirtCapName, dirtCapOffset, dirtInfo.centerInWorld);

    dirtCapInfo.solid = new G4Paraboloid( dirtCapName, capHalfHeight, capTopR, capBottomR);
    
    dirtCapInfo.logical= new G4LogicalVolume( dirtCapInfo.solid, 
                                              dirtMaterial, 
                                              dirtCapName);
    
    dirtCapInfo.physical =  new G4PVPlacement( dirtCapRot, 
                                               dirtCapOffset,
                                               dirtCapInfo.logical,
                                               dirtCapName, 
                                               worldInfo.logical, 
                                               0, 
                                               0,
                                               _config->getBool("g4.doSurfaceCheck",false));

    _visAttributes.push_back(G4VisAttributes(dirtCapVisible, G4Colour::Green()));
    G4VisAttributes& visAtt = _visAttributes.back();
    visAtt.SetForceSolid(dirtCapSolid);
    visAtt.SetForceAuxEdgeVisible(_config->getBool("g4.forceAuxEdgeVisible",false));
    dirtCapInfo.logical->SetVisAttributes(&visAtt);

    addVolInfo( dirtCapInfo );

    return dirtInfo;

  }  // end of Mu2eWorld::constructDirt

  // Construct the DS.  Parent volume is the air inside of the hall.
  // This makes 4 volumes:
  //  0 - a single volume that represents the coils+cryostats in an average way.
  //  1 - DS1, a small piece of DS vacuum that surrounds TS5.
  //  2 - DS2, the upstream part of the DS vacuum, that has a graded field.
  //  3 - DS3, the downstream part of the DS vacuum, that may have a uniform field.
  void Mu2eWorld::constructDS( const VolumeInfo& parent ){

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

    // All Vacuumn volumes fit inside the DS coil+cryostat package.
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

    bool toyDSVisible = _config->getBool("toyDS.visible",true);
    bool toyDSSolid   = _config->getBool("toyDS.solid",true);

    VolumeInfo detSolCoilInfo = nestTubs2( "ToyDSCoil",
                                           detSolCoilParams,
                                           detSolCoilMaterial,
                                           0,
                                           detSolCoilPosition-_hallOriginInMu2e,
                                           parent,
                                           0,
                                           toyDSVisible,
                                           G4Color::Magenta(),
                                           toyDSSolid
                                           );

    // Upstream face of the DS coils+cryo.
    VolumeInfo dsFrontInfo = nestTubs2( "ToyDSFront",
                                        dsFrontParams,
                                        detSolCoilMaterial,
                                        0,
                                        dsFrontPosition-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        toyDSVisible,
                                        G4Color::Blue(),
                                        toyDSSolid
                                        );


    VolumeInfo ds1VacInfo   = nestTubs2( "ToyDS1Vacuum",
                                         ds1VacParams,
                                         vacuumMaterial,
                                         0,
                                         ds1Position-_hallOriginInMu2e,
                                         parent,
                                         0,
                                         toyDSVisible,
                                         G4Colour::Green(),
                                         toyDSSolid
                                         );

    VolumeInfo ds2VacInfo   = nestTubs2( "ToyDS2Vacuum",
                                         ds2VacParams,
                                         vacuumMaterial,
                                         0,
                                         ds2Position-_hallOriginInMu2e,
                                         parent,
                                         0,
                                         toyDSVisible,
                                         G4Colour::Yellow(),
                                         toyDSSolid
                                         );

    VolumeInfo ds3VacInfo = nestTubs2( "ToyDS3Vacuum",
                                       ds3VacParams,
                                       vacuumMaterial,
                                       0,
                                       ds3Position-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       toyDSVisible,
                                       G4Color::Blue(),
                                       toyDSSolid
                                       );

  } // end of Mu2eWorld::constructDS;

  // Construct the TS.  Parent volume is the air inside of the hall.
  void Mu2eWorld::constructTS( const VolumeInfo& parent ){

    // A helper class.
    MaterialFinder materialFinder(*_config);

    // Extract base parameters from config information.
    GeomHandle<Beamline> beamg;
    double solenoidOffset = beamg->solenoidOffset();
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
    vector<double> coll3x, coll3y;
    _config->getVectorDouble("coll3.x",coll3x);
    _config->getVectorDouble("coll3.y",coll3y);

    G4Material* coll1Material  = materialFinder.get("coll1.materialName");
    G4Material* coll3Material  = materialFinder.get("coll3.materialName");
    G4Material* coll5Material  = materialFinder.get("coll5.materialName");

    double coll3Hole           = _config->getDouble("coll3.hole");
    double coll3RotationAngle  = _config->getDouble("coll3.rotationAngle");

    bool collVisible          = _config->getBool("coll.visible",true);
    bool collSolid            = _config->getBool("coll.solid",true);

    // For how all pieces are made from one of two types of material,
    // vacuum or average coils + cryostat material.
    G4Material* vacuumMaterial  = materialFinder.get("toyDS.insideMaterialName");
    G4Material* cryoMaterial    = materialFinder.get("toyDS.materialName");

    GeomHandle<VirtualDetector> vdg;
    double vdHalfLength = vdg->getHalfLength()*CLHEP::mm;

    // Computed quantities.
    double ts5zOffset    = ( rTorus+ts5HalfLength);

    // Build TS1.
    TubsParams ts1VacParams (   0.,  rVac, ts1HalfLength);
    TubsParams ts1CryoParams( rVac, rCryo, ts1HalfLength);

    VolumeInfo ts1VacInfo = nestTubs2( "ToyTS1Vacuum",
                                       ts1VacParams,
                                       vacuumMaterial,
				       beamg->getTS().getTS1().getRotation(),
                                       beamg->getTS().getTS1().getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       toyTSVisible,
                                       G4Color::Red(),
                                       toyTSSolid
                                       );

    VolumeInfo ts1CryoInfo = nestTubs2( "ToyTS1Cryo",
                                        ts1CryoParams,
                                        cryoMaterial,
					beamg->getTS().getTS1().getRotation(),
					beamg->getTS().getTS1().getGlobal()-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        toyTSVisible,
                                        G4Color::Red(),
                                        toyTSSolid
                                        );

    // Place collimator 1

    double coll1Param[7] = { coll1InnerRadius1, rVac, 
			     coll1InnerRadius2, rVac,
			     coll1HalfLength-2*vdHalfLength, 
			     0.0, 360.0*CLHEP::degree };

    VolumeInfo coll1VacInfo = nestCons2( "Coll1",
					 coll1Param,
					 coll1Material,
					 0,
					 beamg->getTS().getColl1().getLocal(),
					 ts1VacInfo,
					 0,
					 collVisible,
					 G4Color::Blue(),
					 collSolid
					 );

    // Build TS2.
    double ts2VacParams[5]  = { 0.0,   rVac, rTorus, 1.5*M_PI, 0.5*M_PI };
    double ts2CryoParams[5] = { rVac, rCryo, rTorus, 1.5*M_PI, 0.5*M_PI };

    // Position in the Mu2e coordintate system.
    G4ThreeVector ts2VacPosition( ts3HalfLength, 0., -rTorus);

    G4RotationMatrix* ts2Rot = new G4RotationMatrix();
    ts2Rot->rotateX(90.0*CLHEP::degree);

    VolumeInfo ts2VacInfo = nestTorus2("ToyTS2Vacuum",
                                       ts2VacParams,
                                       vacuumMaterial,
                                       ts2Rot,
                                       ts2VacPosition-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       toyTSVisible,
                                       G4Color::Yellow(),
                                       toyTSSolid
                                       );

    VolumeInfo ts2CryoInfo = nestTorus2("ToyTS2Cryo",
                                        ts2CryoParams,
                                        cryoMaterial,
                                        ts2Rot,
                                        ts2VacPosition-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        toyTSVisible,
                                        G4Color::Yellow(),
                                        toyTSSolid
                                        );

    // Build TS3.
    TubsParams ts3VacParams (   0.,  rVac, ts3HalfLength);
    TubsParams ts3CryoParams( rVac, rCryo, ts3HalfLength);
    
    VolumeInfo ts3VacInfo = nestTubs2( "ToyTS3Vacuum",
                                       ts3VacParams,
                                       vacuumMaterial,
				       beamg->getTS().getTS3().getRotation(),
				       beamg->getTS().getTS3().getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       toyTSVisible,
                                       G4Color::Green(),
                                       toyTSSolid
                                       );

    VolumeInfo ts3CryoInfo = nestTubs2( "ToyTS3Cryo",
                                        ts3CryoParams,
                                        cryoMaterial,
					beamg->getTS().getTS3().getRotation(),
					beamg->getTS().getTS3().getGlobal()-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        toyTSVisible,
                                        G4Color::Green(),
                                        toyTSSolid
                                        );

    // Place collimator 3

    G4RotationMatrix* coll31Rot = new G4RotationMatrix();
    G4RotationMatrix* coll32Rot = new G4RotationMatrix();
    coll31Rot->rotateZ((180.0+coll3RotationAngle)*CLHEP::degree);
    coll32Rot->rotateZ((180.0+coll3RotationAngle)*CLHEP::degree);
    coll32Rot->rotateY(180.0*CLHEP::degree);

    VolumeInfo coll3Info1 = nestExtrudedSolid2 ( "Coll3_1",
						 coll31HalfLength, coll3x, coll3y,
						 coll3Material,
						 coll31Rot,
						 beamg->getTS().getColl31().getLocal(),
						 ts3VacInfo,
						 0,
						 collVisible,
						 G4Color::Blue(),
						 collSolid
						 );

    VolumeInfo coll3Info2 = nestExtrudedSolid2 ( "Coll3_2",
						 coll31HalfLength, coll3x, coll3y,
						 coll3Material,
						 coll32Rot,
						 beamg->getTS().getColl31().getLocal(),
						 ts3VacInfo,
						 0,
						 collVisible,
						 G4Color::Blue(),
						 collSolid
						 );

    VolumeInfo coll3Info3 = nestExtrudedSolid2 ( "Coll3_3",
						 coll32HalfLength, coll3x, coll3y,
						 coll3Material,
						 coll31Rot,
						 beamg->getTS().getColl32().getLocal(),
						 ts3VacInfo,
						 0,
						 collVisible,
						 G4Color::Blue(),
						 collSolid
						 );

    VolumeInfo coll3Info4 = nestExtrudedSolid2 ( "Coll3_4",
						 coll32HalfLength, coll3x, coll3y,
						 coll3Material,
						 coll32Rot,
						 beamg->getTS().getColl32().getLocal(),
						 ts3VacInfo,
						 0,
						 collVisible,
						 G4Color::Blue(),
						 collSolid
						 );

    double pbarHalfLength     = _config->getDouble("pbar.halfLength");
    G4Material* pbarMaterial  = materialFinder.get("pbar.materialName");
    double pbarParams[5]  = { 0.0,   rVac, pbarHalfLength, 0.0, 360.0*CLHEP::degree };

    VolumeInfo pbarInfo = nestTubs2( "PbarAbs",
				     pbarParams,
				     pbarMaterial,
				     0,
				     G4ThreeVector(0.,0.,0.),
				     ts3VacInfo,
				     0,
				     collVisible,
				     G4Color::Yellow(),
				     collSolid
				     );

    // Build TS4.

    double ts4VacParams[5]  = { 0.0,  rVac, rTorus, 0.5*M_PI, 0.5*M_PI };
    double ts4CryoParams[5] = {rVac, rCryo, rTorus, 0.5*M_PI, 0.5*M_PI };

    // Position in the Mu2e coordintate system.

    G4ThreeVector ts4VacPosition( -ts3HalfLength, 0., rTorus);

    G4RotationMatrix* ts4Rot = new G4RotationMatrix();
    ts4Rot->rotateX(90.0*CLHEP::degree);

    VolumeInfo ts4VacInfo = nestTorus2("ToyTS4Vacuum",
                                       ts4VacParams,
                                       vacuumMaterial,
                                       ts4Rot,
                                       ts4VacPosition-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       toyTSVisible,
                                       G4Color::Yellow(),
                                       toyTSSolid
                                       );

    VolumeInfo ts4CryoInfo = nestTorus2("ToyTS4Cryo",
                                        ts4CryoParams,
                                        cryoMaterial,
                                        ts4Rot,
                                        ts4VacPosition-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        toyTSVisible,
                                        G4Color::Yellow(),
                                        toyTSSolid
                                        );
    
    // Build TS5.

    TubsParams ts5VacParams (   0.,  rVac, ts5HalfLength);
    TubsParams ts5CryoParams( rVac, rCryo, ts5HalfLength);

    VolumeInfo ts5VacInfo = nestTubs2( "ToyTS5Vacuum",
                                       ts5VacParams,
                                       vacuumMaterial,
				       beamg->getTS().getTS5().getRotation(),
				       beamg->getTS().getTS5().getGlobal()-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       toyTSVisible,
                                       G4Color::Red(),
                                       toyTSSolid
                                       );

    VolumeInfo ts5CryoInfo = nestTubs2( "ToyTS5Cryo",
                                        ts5CryoParams,
                                        cryoMaterial,
					beamg->getTS().getTS5().getRotation(),
					beamg->getTS().getTS5().getGlobal()-_hallOriginInMu2e,
                                        parent,
                                        0,
                                        toyTSVisible,
                                        G4Color::Red(),
                                        toyTSSolid
                                        );

    // Place collimator 5

    double coll5Param[5] = { coll5InnerRadius, rVac, 
			     coll5HalfLength-2*vdHalfLength, 0.0, 360.0*CLHEP::degree };

    VolumeInfo coll5VacInfo = nestTubs2( "Coll5",
					 coll5Param,
					 coll5Material,
					 0,
					 beamg->getTS().getColl5().getLocal(),
					 ts5VacInfo,
					 0,
					 collVisible,
					 G4Color::Blue(),
					 collSolid
					 );

  } // end Mu2eWorld::constructTS

  void Mu2eWorld::constructPS( const VolumeInfo& parent ){

    // A helper class.
    MaterialFinder materialFinder(*_config);

    // Extract information from the config file.
    
    GeomHandle<Beamline> beamg;
    double solenoidOffset = beamg->solenoidOffset();
    double rTorus         = beamg->getTS().torusRadius();
    double ts1HalfLength  = beamg->getTS().getTS1().getHalfLength();

    double ps1HalfLength     = _config->getDouble("toyPS1.vacHalfLength");

    // Build the barrel of the cryostat.
    TubsParams psCryoParams( _config->getDouble("toyPS.rIn"),
                             _config->getDouble("toyPS.rOut"),
                             _config->getDouble("toyPS.CryoHalfLength"));
    G4Material* psCryoMaterial = materialFinder.get("toyPS.materialName");

    // In the Mu2e coordinate system.
    double psCryoZ0 = -rTorus + -2.*ts1HalfLength - psCryoParams.zHalfLength;
    G4ThreeVector psCryoPosition( solenoidOffset, 0., psCryoZ0 );
    
    bool toyPSVisible = _config->getBool("toyPS.visible",true);
    bool toyPSSolid   = _config->getBool("toyPS.solid",true);

    
    // Toy model of the PS coils + cryostat. It needs real structure.
    VolumeInfo psCryoInfo = nestTubs2( "PSCryo",
                                       psCryoParams,
                                       psCryoMaterial,
                                       0,
                                       psCryoPosition-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       toyPSVisible,
                                       G4Color::Cyan(),
                                       toyPSSolid
                                       );

    // Build the main PS vacuum body.
    TubsParams ps1VacParams( 0.,
                             _config->getDouble("toyPS.rIn"),
                             ps1HalfLength);
    G4Material* vacuumMaterial  = materialFinder.get("toyDS.insideMaterialName");

    // Position in the Mu2e coordinate system.
    double ps1Z0     = -rTorus + -2.*ts1HalfLength - ps1HalfLength;
    G4ThreeVector ps1Position( solenoidOffset, 0., ps1Z0);    
    
    VolumeInfo ps1VacInfo   = nestTubs2( "PS1Vacuum",
                                         ps1VacParams,
                                         vacuumMaterial,
                                         0,
                                         ps1Position-_hallOriginInMu2e,
                                         parent,
                                         0,
                                         toyPSVisible,
                                         G4Colour::Green(),
                                         toyPSSolid
                                         );
    
    // Build the production target.
    TubsParams prodTargetParams( 0.,
                                 _config->getDouble("targetPS_rOut"),
                                 _config->getDouble("targetPS_halfLength"));
    G4Material* prodTargetMaterial = materialFinder.get("targetPS_materialName");
    
    // Position in the Mu2e coordinate system.
    CLHEP::Hep3Vector prodTargetPosition = _config->getHep3Vector("productionTarget.position");

    // Rotation of production target.
    double targetPS_rotX = _config->getDouble("targetPS_rotX" )*CLHEP::degree;
    double targetPS_rotY = _config->getDouble("targetPS_rotY" )*CLHEP::degree;

    // G4 takes ownership of this G4RotationMatrix object.
    // Passive rotation. See Mu2e-doc-938.
    G4RotationMatrix* prodTargetRotation = new G4RotationMatrix;
    prodTargetRotation->rotateY( -targetPS_rotY);
    prodTargetRotation->rotateX( -targetPS_rotX);

    bool prodTargetVisible = _config->getBool("targetPS.visible",true);
    bool prodTargetSolid   = _config->getBool("targetPS.solid",true);

    VolumeInfo prodTargetInfo   = nestTubs2( "ProductionTarget",
                                             prodTargetParams,
                                             prodTargetMaterial,
                                             prodTargetRotation,
                                             prodTargetPosition-ps1Position,
                                             ps1VacInfo,
                                             0,
                                             prodTargetVisible,
                                             G4Colour::Magenta(),
                                             prodTargetSolid
                                             );

    
    // Set the parameters of the transformation from the PrimaryProtonGun
    // coordinates to G4 coordinates.  This needs an active sense rotation,
    // the opposite of what G4 needed.
    _primaryProtonGunRotation = prodTargetRotation->inverse();

    G4ThreeVector prodTargetFaceLocal(0.,0.,prodTargetParams.zHalfLength);
    _primaryProtonGunOrigin = prodTargetPosition + _mu2eOrigin + _primaryProtonGunRotation*prodTargetFaceLocal;
    
  } // end Mu2eWorld::constructPS


  // Choose the selected tracker and build it.
  VolumeInfo Mu2eWorld::constructTracker(){

    // The tracker is built inside this volume.
    VolumeInfo detSolDownstreamVacInfo = locateVolInfo("ToyDS3Vacuum");

    // z Position of the center of the DS solenoid parts, given in the Mu2e coordinate system.
    double z0DSdown = detSolDownstreamVacInfo.centerInWorld.z()+_hallOriginInMu2e.z();

    // Construct one of the trackers.
    VolumeInfo trackerInfo;
    if( _config->getBool("hasLTracker",false) ){
      int ver = _config->getInt("LTrackerVersion",3);
      //cout << "LTracker version: " << ver << "\n";
      if ( ver == 3 ){
        trackerInfo = constructLTrackerv3( detSolDownstreamVacInfo.logical, z0DSdown, *_config );
      }
    } else if ( _config->getBool("hasITracker",false) ) {
      trackerInfo = ITrackerBuilder::constructTracker( detSolDownstreamVacInfo.logical, z0DSdown );
    } else if ( _config->getBool("hasTTracker",false) ) {
      int ver = _config->getInt("TTrackerVersion",3);
      if ( ver == 3 ){
        trackerInfo = constructTTrackerv3( detSolDownstreamVacInfo.logical, z0DSdown, *_config );
      }        
    } else {
      trackerInfo = constructDummyTracker( detSolDownstreamVacInfo.logical, z0DSdown, *_config );
    }

    return trackerInfo;

  } // end Mu2eWorld::constructTracker

  // Either build the stopping target or a placeholder.
  // This should be call constrcutStoppingTarget but the name is already taken.
  VolumeInfo Mu2eWorld::constructTarget(){


    // The target is built inside this volume.
    VolumeInfo detSolUpstreamVacInfo   = locateVolInfo("ToyDS2Vacuum");

    // z Position of the center of the DS solenoid parts, given in the Mu2e coordinate system.
    double z0DSup   = detSolUpstreamVacInfo.centerInWorld.z()+_hallOriginInMu2e.z();

    // Buid the stopping target
    VolumeInfo targetInfo;
    if( _config->getBool("hasTarget",false) ){

      targetInfo = constructStoppingTarget( detSolUpstreamVacInfo.logical, 
                                            z0DSup,
                                            *_config );

    } else {

      targetInfo = constructDummyStoppingTarget( detSolUpstreamVacInfo.logical, 
                                                 z0DSup,
                                                 *_config );
    } //hasTarget

    return targetInfo;
  } // end Mu2eWorld::constructTarget

  //construct proton absorber
  void Mu2eWorld::constructProtonAbs(){
    
    VolumeInfo parent1 = locateVolInfo("ToyDS2Vacuum");
    VolumeInfo parent2 = locateVolInfo("ToyDS3Vacuum");
    double pabs1rIn0   = _config->getDouble("protonabsorber.InRadius0");
    double pabs1rOut0  = _config->getDouble("protonabsorber.OutRadius0");
    double pabs2rIn1   = _config->getDouble("protonabsorber.InRadius1");
    double pabs2rOut1  = _config->getDouble("protonabsorber.OutRadius1");
    double zLen   = _config->getDouble("protonabsorber.halfLength");
    double thick  = _config->getDouble("protonabsorber.thickness");
    double trkMinr= _config->getDouble("ttracker.envelopeInnerRadius");
    pabs1rIn0 = pabs1rOut0 - thick;
    pabs2rOut1 = trkMinr;  pabs2rIn1 = trkMinr - thick;
    
    MaterialFinder materialFinder(*_config);
    G4Material* pabsMaterial = materialFinder.get("protonabsorber.materialName");
    
    double z0DSup   = parent1.centerInWorld.z()+_hallOriginInMu2e.z();
    vector<double> targetRadius;  _config->getVectorDouble("target.radii", targetRadius);
    double numoftf = (targetRadius.size()-1)/2.0;
    double foilwid=_config->getDouble("target.deltaZ"); double taglen =(foilwid*numoftf) + 5.0;
    double z0valt =_config->getDouble("target.z0");     double tagoff =z0valt - z0DSup + 12000.0;
    double targetEnd = tagoff + taglen;
    double ds2halflen = _config->getDouble("toyDS2.halfLength");
    double pabs1len = ds2halflen - targetEnd;
    G4ThreeVector  pabs1Offset(0.0, 0.0, (pabs1len/2.0) + targetEnd);
    
    double pabs1rOut1 = ((pabs2rOut1 - pabs1rOut0)*(pabs1len/(2.0*zLen))) + pabs1rOut0;
    double pabs1rIn1  = pabs1rOut1 - thick;
    double ds3halflen = _config->getDouble("toyDS3.halfLength");
    double pabs2len  = (2.0*zLen) - pabs1len;
    double pabs2zoff = (pabs2len/2.0) - ds3halflen;
    G4ThreeVector  pabs2Offset(0.0, 0.0, pabs2zoff);
    
    // proton absorber in DS2
    double pabs1Param[7] = { pabs1rIn0, pabs1rOut0, pabs1rIn1, pabs1rOut1, pabs1len/2.0, 0.0, 360.0*CLHEP::degree };
    bool pabsVisible = _config->getBool("protonabsorber.visible",true);
    bool pabsSolid   = _config->getBool("protonabsorber.solid",true);
    
    if( _config->getBool("hasProtonAbsorber", true) ){
      
      edm::LogInfo log("GEOM");
      log << "Constructing Proton Absorber -- \n";
      log << "Proton Abs Offset in DS2:  " << pabs1Offset <<"\n";
      log << "rIn,  rOut (-z): "<< pabs1rIn0 <<"  "<< pabs1rOut0<<"  ";
      log << "rIn,  rOut (+z): "<< pabs1rIn1 <<"  "<< pabs1rOut1<<"  ";
      log << "halflength: "<< pabs1len/2.0 <<"\n";
      log << "Proton Abs Offset in DS3:  " << pabs2Offset <<"\n";
      log << "rIn,  rOut (-z): "<< pabs1rIn1 <<"  "<< pabs1rOut1<<"  ";
      log << "rIn,  rOut (+z): "<< pabs2rIn1 <<"  "<< pabs2rOut1<<"  ";
      log << "halflength: "<< pabs2len/2.0 <<"\n";
      
      VolumeInfo protonabs1Info = nestCons2( "protonabs1",
					     pabs1Param,
					     pabsMaterial,
					     0,
					     pabs1Offset,
					     parent1,
					     0,
					     pabsVisible,
					     G4Color::White(),
					     pabsSolid
					     );
      
      // proton absorber in DS3
      double pabs2Param[7] = { pabs1rIn1, pabs1rOut1, pabs2rIn1, pabs2rOut1, pabs2len/2.0, 0.0, 360.0*CLHEP::degree };
      
      VolumeInfo protonabs2Info = nestCons2( "protonabs2",
					     pabs2Param,
					     pabsMaterial,
					     0,
					     pabs2Offset,
					     parent2,
					     0,
					     pabsVisible,
					     G4Color::Magenta(),
					     pabsSolid
					     );
    }
  } // end of Mu2eWorld::constructProtonAbs;

  // Construct the magnetic field managers and attach them to
  // the relevant volumes.
  void Mu2eWorld::constructBFieldAndManagers(){

    // Figure out which magnetic field managers are needed.
    int dsFieldForm    = _config->getInt("detSolFieldForm", dsModelUniform); 

    // Decide on the G4 Stepper

    bool needDSFull    = (dsFieldForm == dsModelFull  || dsFieldForm == dsModelSplit );
    bool needDSUniform = (dsFieldForm == dsModelSplit || dsFieldForm == dsModelUniform );

    string stepper = _config->getString("g4.stepper","G4SimpleRunge");

    // Create field manager for the full DS field.
    if ( needDSFull ){
      if ( stepper  == "G4ClassicalRK4" ) {
        _dsFull = FieldMgr::forMappedField<G4ClassicalRK4>( "DS", _mu2eOrigin );
      } else if ( stepper  == "G4ImplicitEuler" ) {
        _dsFull = FieldMgr::forMappedField<G4ImplicitEuler>( "DS", _mu2eOrigin );
      } else if ( stepper  == "G4ExplicitEuler" ) {
        _dsFull = FieldMgr::forMappedField<G4ExplicitEuler>( "DS", _mu2eOrigin );
      } else if ( stepper  == "G4SimpleHeum" ) {
        _dsFull = FieldMgr::forMappedField<G4SimpleHeum>( "DS", _mu2eOrigin );
      } else if ( stepper  == "G4HelixImplicitEuler" ) {
        _dsFull = FieldMgr::forMappedField<G4HelixImplicitEuler>( "DS", _mu2eOrigin );
      } else if ( stepper  == "G4HelixSimpleRunge" ) {
        _dsFull = FieldMgr::forMappedField<G4HelixSimpleRunge>( "DS", _mu2eOrigin );
      } else {
        _dsFull = FieldMgr::forMappedField<G4SimpleRunge>( "DS", _mu2eOrigin );
      } 
    }

    // Create field manager for the uniform DS field.
    if ( needDSUniform){
      // Handle to the BField manager.
      GeomHandle<BFieldManager> bfMgr;
      _dsUniform = FieldMgr::forUniformField( bfMgr->getDSUniformValue(), _mu2eOrigin );
    }

    // Create field managers for the PS and TS.

    if ( stepper  == "G4ClassicalRK4" ) {
      _psFull = FieldMgr::forMappedField<G4ClassicalRK4>( "PS", _mu2eOrigin );
    } else if ( stepper  == "G4ImplicitEuler" ) {
      _psFull = FieldMgr::forMappedField<G4ImplicitEuler>( "PS", _mu2eOrigin );
    } else if ( stepper  == "G4ExplicitEuler" ) {
      _psFull = FieldMgr::forMappedField<G4ExplicitEuler>( "PS", _mu2eOrigin );
    } else if ( stepper  == "G4SimpleHeum" ) {
      _psFull = FieldMgr::forMappedField<G4SimpleHeum>( "PS", _mu2eOrigin );
    } else if ( stepper  == "G4HelixImplicitEuler" ) {
      _psFull = FieldMgr::forMappedField<G4HelixImplicitEuler>( "PS", _mu2eOrigin );
    } else if ( stepper  == "G4HelixSimpleRunge" ) {
      _psFull = FieldMgr::forMappedField<G4HelixSimpleRunge>( "PS", _mu2eOrigin );
    } else {
      _psFull = FieldMgr::forMappedField<G4SimpleRunge>( "PS", _mu2eOrigin );
    } 


    if ( stepper  == "G4ClassicalRK4" ) {
      _tsFull = FieldMgr::forMappedField<G4ClassicalRK4>( "TS", _mu2eOrigin );
    } else if ( stepper  == "G4ImplicitEuler" ) {
      _tsFull = FieldMgr::forMappedField<G4ImplicitEuler>( "TS", _mu2eOrigin );
    } else if ( stepper  == "G4ExplicitEuler" ) {
      _tsFull = FieldMgr::forMappedField<G4ExplicitEuler>( "TS", _mu2eOrigin );
    } else if ( stepper  == "G4SimpleHeum" ) {
      _tsFull = FieldMgr::forMappedField<G4SimpleHeum>( "TS", _mu2eOrigin );
    } else if ( stepper  == "G4HelixImplicitEuler" ) {
      _tsFull = FieldMgr::forMappedField<G4HelixImplicitEuler>( "TS", _mu2eOrigin );
    } else if ( stepper  == "G4HelixSimpleRunge" ) {
      _tsFull = FieldMgr::forMappedField<G4HelixSimpleRunge>( "TS", _mu2eOrigin );
    } else {
      _tsFull = FieldMgr::forMappedField<G4SimpleRunge>( "TS", _mu2eOrigin );
    }

    // Get pointers to logical volumes.
    G4LogicalVolume* ds2Vacuum = locateVolInfo("ToyDS2Vacuum").logical;
    G4LogicalVolume* ds3Vacuum = locateVolInfo("ToyDS3Vacuum").logical;

    vector<G4LogicalVolume*> psVacua;
    psVacua.push_back( locateVolInfo("PS1Vacuum").logical );

    vector<G4LogicalVolume*> tsVacua;
    tsVacua.push_back( locateVolInfo("ToyTS1Vacuum").logical );
    tsVacua.push_back( locateVolInfo("ToyTS2Vacuum").logical );
    tsVacua.push_back( locateVolInfo("ToyTS3Vacuum").logical );
    tsVacua.push_back( locateVolInfo("ToyTS4Vacuum").logical );
    tsVacua.push_back( locateVolInfo("ToyTS5Vacuum").logical );

    // Attach field managers to the appropriate logical volumes.
    if (dsFieldForm == dsModelFull  ){
      ds2Vacuum->SetFieldManager( _dsFull->manager(), true);
      ds3Vacuum->SetFieldManager( _dsFull->manager(), true);
    } else if ( dsFieldForm == dsModelSplit ){
      ds2Vacuum->SetFieldManager( _dsFull->manager(), true);
      ds3Vacuum->SetFieldManager( _dsUniform->manager(), true);
    } else {
      ds2Vacuum->SetFieldManager( _dsUniform->manager(), true);
      ds3Vacuum->SetFieldManager( _dsUniform->manager(), true);
    }

    for ( vector<G4LogicalVolume*>::iterator i=psVacua.begin();
          i!=psVacua.end(); ++i ){
      (**i).SetFieldManager( _psFull->manager(), true);
    }
    
    for ( vector<G4LogicalVolume*>::iterator i=tsVacua.begin();
          i!=tsVacua.end(); ++i ){
      (**i).SetFieldManager( _tsFull->manager(), true);
    }

    // Adjust properties of the integrators to control accuracy vs time.
    G4double singleValue         = 0.5e-01*CLHEP::mm;
    G4double newUpstreamDeltaI   = singleValue;
    G4double newDownstreamDeltaI = singleValue;
    G4double deltaOneStep        = singleValue;
    G4double deltaChord          = singleValue;

    if ( _dsFull.get() != 0 ){
      _dsFull->manager()->SetDeltaOneStep(deltaOneStep);
      _dsFull->manager()->SetDeltaIntersection(newUpstreamDeltaI);
      _dsFull->chordFinder()->SetDeltaChord(deltaChord);
    }

    if ( _dsUniform.get() != 0 ){
      G4double deltaIntersection = 0.00001*CLHEP::mm;
      _dsUniform->manager()->SetDeltaIntersection(deltaIntersection);
    }

    _tsFull->manager()->SetDeltaOneStep(deltaOneStep);
    _tsFull->manager()->SetDeltaIntersection(newUpstreamDeltaI);
    _tsFull->chordFinder()->SetDeltaChord(deltaChord);

    _psFull->manager()->SetDeltaOneStep(deltaOneStep);
    _psFull->manager()->SetDeltaIntersection(newUpstreamDeltaI);
    _psFull->chordFinder()->SetDeltaChord(deltaChord);


  } // end Mu2eWorld::constructBFieldAndManagers

  // Construct the magnetic field managers and attach them to
  // the relevant volumes - alternative version.
  void Mu2eWorld::constructBFieldAndManagers2(){

    // Figure out which magnetic field managers are needed.
    int dsFieldForm    = _config->getInt("detSolFieldForm", dsModelUniform); 

    // Decide on the G4 Stepper

    bool needDSUniform = (dsFieldForm == dsModelSplit || dsFieldForm == dsModelUniform );

    string stepper = _config->getString("g4.stepper","G4SimpleRunge");

    // Create field manager for the uniform DS field.
    if ( needDSUniform){
      // Handle to the BField manager.
      GeomHandle<BFieldManager> bfMgr;
      _dsUniform = FieldMgr::forUniformField( bfMgr->getDSUniformValue(), _mu2eOrigin );
    }

    // Create field managers for the PS and TS.

    if ( stepper  == "G4ClassicalRK4" ) {
      _hallFull = FieldMgr::forMappedField<G4ClassicalRK4>( "", _mu2eOrigin );
    } else if ( stepper  == "G4ImplicitEuler" ) {
      _hallFull = FieldMgr::forMappedField<G4ImplicitEuler>( "", _mu2eOrigin );
    } else if ( stepper  == "G4ExplicitEuler" ) {
      _hallFull = FieldMgr::forMappedField<G4ExplicitEuler>( "", _mu2eOrigin );
    } else if ( stepper  == "G4SimpleHeum" ) {
      _hallFull = FieldMgr::forMappedField<G4SimpleHeum>( "", _mu2eOrigin );
    } else if ( stepper  == "G4HelixImplicitEuler" ) {
      _hallFull = FieldMgr::forMappedField<G4HelixImplicitEuler>( "", _mu2eOrigin );
    } else if ( stepper  == "G4HelixSimpleRunge" ) {
      _hallFull = FieldMgr::forMappedField<G4HelixSimpleRunge>( "", _mu2eOrigin );
    } else {
      _hallFull = FieldMgr::forMappedField<G4SimpleRunge>( "", _mu2eOrigin );
    } 

    // Get pointers to logical volumes.
    G4LogicalVolume* ds2Vacuum = locateVolInfo("ToyDS2Vacuum").logical;
    G4LogicalVolume* ds3Vacuum = locateVolInfo("ToyDS3Vacuum").logical;

    vector<G4LogicalVolume*> psVacua;
    psVacua.push_back( locateVolInfo("PS1Vacuum").logical );

    vector<G4LogicalVolume*> tsVacua;
    tsVacua.push_back( locateVolInfo("ToyTS1Vacuum").logical );
    tsVacua.push_back( locateVolInfo("ToyTS2Vacuum").logical );
    tsVacua.push_back( locateVolInfo("ToyTS3Vacuum").logical );
    tsVacua.push_back( locateVolInfo("ToyTS4Vacuum").logical );
    tsVacua.push_back( locateVolInfo("ToyTS5Vacuum").logical );

    // Attach field managers to the appropriate logical volumes.
    if (dsFieldForm == dsModelFull  ){
      ds2Vacuum->SetFieldManager( _hallFull->manager(), true);
      ds3Vacuum->SetFieldManager( _hallFull->manager(), true);
    } else if ( dsFieldForm == dsModelSplit ){
      ds2Vacuum->SetFieldManager( _hallFull->manager(), true);
      ds3Vacuum->SetFieldManager( _dsUniform->manager(), true);
    } else {
      ds2Vacuum->SetFieldManager( _dsUniform->manager(), true);
      ds3Vacuum->SetFieldManager( _dsUniform->manager(), true);
    }

    for ( vector<G4LogicalVolume*>::iterator i=psVacua.begin();
          i!=psVacua.end(); ++i ){
      (**i).SetFieldManager( _hallFull->manager(), true);
    }
    
    for ( vector<G4LogicalVolume*>::iterator i=tsVacua.begin();
          i!=tsVacua.end(); ++i ){
      (**i).SetFieldManager( _hallFull->manager(), true);
    }

    // Adjust properties of the integrators to control accuracy vs time.
    G4double singleValue         = 0.5e-01*CLHEP::mm;
    G4double newUpstreamDeltaI   = singleValue;
    G4double newDownstreamDeltaI = singleValue;
    G4double deltaOneStep        = singleValue;
    G4double deltaChord          = singleValue;

    if ( _dsUniform.get() != 0 ){
      G4double deltaIntersection = 0.00001*CLHEP::mm;
      _dsUniform->manager()->SetDeltaIntersection(deltaIntersection);
    }

    _hallFull->manager()->SetDeltaOneStep(deltaOneStep);
    _hallFull->manager()->SetDeltaIntersection(newUpstreamDeltaI);
    _hallFull->chordFinder()->SetDeltaChord(deltaChord);

  } // end Mu2eWorld::constructBFieldAndManagers

  // Adding a step limiter is a two step process.
  // 1) In the physics list constructor add a G4StepLimiter to the list of discrete
  //    physics processes attached to each particle species of interest.
  //
  // 2) In this code, create a G4UserLimits object and attach it to the logical 
  //    volumes of interest.
  // The net result is specifying a step limiter for pairs of (logical volume, particle species).
  //
  void Mu2eWorld::constructStepLimiters(){

    // Maximum step length, in mm.
    double maxStep = _config->getDouble("bfield.maxStep", 20.);

    G4LogicalVolume* ds2Vacuum      = locateVolInfo("ToyDS2Vacuum").logical;
    G4LogicalVolume* ds3Vacuum      = locateVolInfo("ToyDS3Vacuum").logical;
    G4LogicalVolume* tracker        = locateVolInfo("TrackerMother").logical;
    G4LogicalVolume* stoppingtarget = locateVolInfo("StoppingTargetMother").logical;

    vector<G4LogicalVolume*> psVacua;
    psVacua.push_back( locateVolInfo("PS1Vacuum").logical );

    vector<G4LogicalVolume*> tsVacua;
    tsVacua.push_back( locateVolInfo("ToyTS1Vacuum").logical );
    tsVacua.push_back( locateVolInfo("ToyTS2Vacuum").logical );
    tsVacua.push_back( locateVolInfo("ToyTS3Vacuum").logical );
    tsVacua.push_back( locateVolInfo("ToyTS4Vacuum").logical );
    tsVacua.push_back( locateVolInfo("ToyTS5Vacuum").logical );

    // We may make separate G4UserLimits objects per logical volume but we choose not to.
    _stepLimits.push_back( G4UserLimits(maxStep) );
    G4UserLimits* stepLimit = &(_stepLimits.back());

    ds2Vacuum->SetUserLimits( stepLimit );
    ds3Vacuum->SetUserLimits( stepLimit );

    tracker->SetUserLimits( stepLimit );
    stoppingtarget->SetUserLimits( stepLimit );

    for ( vector<G4LogicalVolume*>::iterator i=psVacua.begin();
          i!=psVacua.end(); ++i ){
      (**i).SetUserLimits( stepLimit);
    }

    for ( vector<G4LogicalVolume*>::iterator i=tsVacua.begin();
          i!=tsVacua.end(); ++i ){
      (**i).SetUserLimits( stepLimit);
    }


  } // end Mu2eWorld::constructStepLimiters(){


  // Construct calorimeter if needed.
  void Mu2eWorld::constructCal(){

    if ( ! _config->getBool("hasCalorimeter",false) ) return;

    VolumeInfo detSolDownstreamVacInfo = locateVolInfo("ToyDS3Vacuum");

    double z0DSdown = detSolDownstreamVacInfo.centerInWorld.z()+_hallOriginInMu2e.z();

    constructCalorimeter( detSolDownstreamVacInfo.logical,
			  -z0DSdown,
			  *_config );
  }


  // A place holder for now.
  void Mu2eWorld::constructMagnetYoke(){
  } // end Mu2eWorld::constructMagnetYoke

  // A place holder for now.
  void Mu2eWorld::constructCRV(){
  } // end Mu2eWorld::constructCRV

  void Mu2eWorld::constructSteel( const VolumeInfo& parent ){

    MaterialFinder materialFinder(*_config);

    // Extract information from the config file
    double HallSteelHalfThick   = _config->getDouble("fluxcrv.HallSteelHalfThick");
    double HallSteelHalfLenXY = _config->getDouble("fluxcrv.HallSteelHalfLengthXY");
    double HallSteelHalfLenZ = _config->getDouble("fluxcrv.HallSteelHalfLengthZ");
    G4Material* HallSteelShieldMaterial = materialFinder.get("fluxcrv.HallSteelMaterialName");

    GeomHandle<Beamline> beamg;
    double rCryo = beamg->getTS().outerRadius();

    G4Material* FrontHoleMaterial = materialFinder.get("toyDS.insideMaterialName");

    // Compute dimensions of 5 sides in Mu2e coordinates
    double HallSteelTopHalfX   = HallSteelHalfLenXY + HallSteelHalfThick;
    double HallSteelTopHalfY   = HallSteelHalfThick;
    double HallSteelTopHalfZ   = HallSteelHalfLenZ;
    double HallSteelSideHalfX  = HallSteelHalfThick;
    double HallSteelSideHalfY  = HallSteelHalfLenXY - HallSteelHalfThick;
    double HallSteelSideHalfZ  = HallSteelHalfLenZ; 
    double HallSteelFrontHalfX = HallSteelHalfLenXY - HallSteelHalfThick;
    double HallSteelFrontHalfY = HallSteelHalfLenXY - HallSteelHalfThick;
    double HallSteelFrontHalfZ = HallSteelHalfThick;

    double HallSteelTopDims[3] ={
      HallSteelTopHalfX,
      HallSteelTopHalfY,
      HallSteelTopHalfZ
    };

    double HallSteelSideDims[3] ={
      HallSteelSideHalfX,
      HallSteelSideHalfY,
      HallSteelSideHalfZ
    };

    double HallSteelFrontDims[3] ={
      HallSteelFrontHalfX,
      HallSteelFrontHalfY,           
      HallSteelFrontHalfZ                        
    };
    TubsParams FrontHoleDims(0.,
                             _config->getDouble("fluxcrv.HallSteelHoleRadius"),
                             HallSteelHalfThick
                             );

    // Get positions of each side. Assuming view from target foils 
    double dsCoilZ0          = _config->getDouble("toyDS.z0");

    double solenoidOffset = beamg->solenoidOffset();

    //G4ThreeVector detSolCoilPosition(-solenoidOffset, 0., -dsCoilZ0);
    G4ThreeVector detSolCoilPosition(+solenoidOffset, 0., -dsCoilZ0);

    std::vector<double> HallSteelOffsetSTDV;
    _config->getVectorDouble("fluxcrv.HallSteelOffset", HallSteelOffsetSTDV, 3);
    G4ThreeVector HallSteelOffset(HallSteelOffsetSTDV[0],
                                  HallSteelOffsetSTDV[1],
                                  HallSteelOffsetSTDV[2]);

    _hallOriginInMu2e = parent.centerInWorld - _mu2eOrigin + detSolCoilPosition + HallSteelOffset;

    G4ThreeVector TopShield   (0.,      HallSteelSideHalfY + HallSteelHalfThick, 0.);
    G4ThreeVector BottomShield(0.,    -(HallSteelSideHalfY + HallSteelHalfThick), 0.);
    G4ThreeVector LeftShield  (         HallSteelSideHalfY + HallSteelHalfThick,0., 0.);
    G4ThreeVector RightShield (       -(HallSteelSideHalfY + HallSteelHalfThick),0., 0.);
    G4ThreeVector BackShield  (0., 0.,  HallSteelSideHalfZ - HallSteelHalfThick);
    G4ThreeVector FrontShield (0., 0.,-(HallSteelSideHalfZ - HallSteelHalfThick));

    //Hole in front shield for TS is centered in the shield
    G4ThreeVector FrontHole(0.,0.,0.);

    bool hallSteelVisible = _config->getBool("fluxcrv.visible",true);
    bool hallSteelSolid   = _config->getBool("fluxcrv.solid",false);

    // Place Boxes

    VolumeInfo TopInfo = nestBox ("HallSteelTopShield",
                                  HallSteelTopDims,
                                  HallSteelShieldMaterial,
                                  0,
                                  TopShield -_hallOriginInMu2e,
                                  parent,
                                  0,
                                  hallSteelVisible,
                                  G4Colour::Green(),
                                  hallSteelSolid);

    VolumeInfo BottomInfo = nestBox ("HallSteelBottomShield",
                                     HallSteelTopDims, 
                                     HallSteelShieldMaterial,
                                     0,
                                     BottomShield -_hallOriginInMu2e,
                                     parent,
                                     0, 
                                     hallSteelVisible,
                                     G4Colour::Green(), 
                                     hallSteelSolid);

    VolumeInfo LeftInfo = nestBox ("HallSteelLeftShield",
                                   HallSteelSideDims,
                                   HallSteelShieldMaterial,
                                   0, 
                                   LeftShield -_hallOriginInMu2e,
                                   parent,
                                   0, 
                                   hallSteelVisible,
                                   G4Colour::Green(),
                                   hallSteelSolid); 

    VolumeInfo RightInfo = nestBox ("HallSteelRightShield",
                                    HallSteelSideDims,
                                    HallSteelShieldMaterial,
                                    0, 
                                    RightShield -_hallOriginInMu2e,
                                    parent,
                                    0, 
                                    hallSteelVisible,
                                    G4Colour::Green(),
                                    hallSteelSolid); 

    VolumeInfo BackInfo = nestBox ("HallSteelBackShield",
                                   HallSteelFrontDims,
                                   HallSteelShieldMaterial,
                                   0, 
                                   BackShield -_hallOriginInMu2e,
                                   parent,
                                   0, 
                                   hallSteelVisible,
                                   G4Colour::Green(),
                                   hallSteelSolid); 

    // constructing "hollow" front shield

    G4Box* CRVBox = new G4Box("CRVFrontShieldBox", 
                              HallSteelFrontDims[0], 
                              HallSteelFrontDims[1], 
                              HallSteelFrontDims[2]);

    G4Tubs* HallSteelFrontHoleTubs = new G4Tubs("HallSteelFrontHoleTubs", 
                                                FrontHoleDims.innerRadius,
                                                FrontHoleDims.outerRadius,
                                                FrontHoleDims.zHalfLength,
                                                FrontHoleDims.phi0, 
                                                FrontHoleDims.phiMax);

    string const FrontShieldName = "CRVFrontShield";
    VolumeInfo FrontShieldInfo;

    FrontShieldInfo.solid = 
      new G4SubtractionSolid(FrontShieldName, CRVBox, HallSteelFrontHoleTubs);

    FrontShieldInfo.logical = new G4LogicalVolume( FrontShieldInfo.solid, 
                                                   HallSteelShieldMaterial, 
                                                   FrontShieldName); 

    FrontShieldInfo.physical = new G4PVPlacement( 0,
                                                  FrontShield -_hallOriginInMu2e, 
                                                  FrontShieldInfo.logical,
                                                  FrontShieldName,
                                                  parent.logical,
                                                  0,
                                                  0,
                                                  _config->getBool("g4.doSurfaceCheck",false));
    
    if ( hallSteelVisible ){
      // We need to manage the lifetime of the G4VisAttributes object.
      _visAttributes.push_back(G4VisAttributes(true, G4Colour::Green()));
      G4VisAttributes& visAtt = _visAttributes.back();
      visAtt.SetForceAuxEdgeVisible(_config->getBool("g4.forceAuxEdgeVisible",false));
      visAtt.SetForceSolid(hallSteelSolid);
      FrontShieldInfo.logical->SetVisAttributes(&visAtt);
    } 
    else{
      FrontShieldInfo.logical->SetVisAttributes(G4VisAttributes::Invisible);
    }

  } // end Mu2eWorld::constructSteel

  // Construct the virtual detectors

  void Mu2eWorld::constructVD( ){

    // Place virtual detectors

    G4SDManager* SDman   = G4SDManager::GetSDMpointer();
    VirtualDetectorSD* vdSD     = new VirtualDetectorSD("VirtualDetector", *_config);
    SDman->AddNewDetector(vdSD);

    bool   vdVisible    = _config->getBool("vd.visible",true);
    bool   vdSolid      = _config->getBool("vd.solid",true);
    
    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()<=0 ) return;

    GeomHandle<Beamline> beamg;
    double rVac           = beamg->getTS().innerRadius();

    double vdHalfLength = vdg->getHalfLength();

    MaterialFinder materialFinder(*_config);
    G4Material* vacuumMaterial     = materialFinder.get("toyDS.insideMaterialName");

    TubsParams vdParams(0,rVac,vdHalfLength);

    // VD 1 and 2 are placed inside TS1

    for( int id=1; id<=2; ++id) if( vdg->exist(id) ) {
      VolumeInfo parent = locateVolInfo("ToyTS1Vacuum");
      ostringstream name;
      name << "VirtualDetector" << id;
      VolumeInfo vd = nestTubs2( name.str(), vdParams, vacuumMaterial, 0,
				 vdg->getLocal(id),
				 parent,
				 id, vdVisible, G4Color::Red(), vdSolid );
      vd.logical->SetSensitiveDetector(vdSD);
    }

    // VD 3-6 are placed inside TS3

    for( int id=3; id<=6; ++id) if( vdg->exist(id) ) {
      VolumeInfo parent = locateVolInfo("ToyTS3Vacuum");
      ostringstream name;
      name << "VirtualDetector" << id;
      VolumeInfo vd = nestTubs2( name.str(), vdParams, vacuumMaterial, 0,
				 vdg->getLocal(id),
				 parent,
				 id, vdVisible, G4Color::Red(), vdSolid );
      vd.logical->SetSensitiveDetector(vdSD);
    }

    // VD 7-8 are placed inside TS3

    for( int id=7; id<=8; ++id) if( vdg->exist(id) ) {
      VolumeInfo parent = locateVolInfo("ToyTS5Vacuum");
      ostringstream name;
      name << "VirtualDetector" << id;
      VolumeInfo vd = nestTubs2( name.str(), vdParams, vacuumMaterial, 0,
				 vdg->getLocal(id),
				 parent,
				 id, vdVisible, G4Color::Red(), vdSolid );
      vd.logical->SetSensitiveDetector(vdSD);
    }

  } //constructVD

  // Place a G4Box inside a logical volume.
  VolumeInfo Mu2eWorld::nestBox ( string const& name,
                                  double const halfDim[3],
                                  G4Material* material,
                                  G4RotationMatrix* rot,
                                  G4ThreeVector const& offset,
                                  const VolumeInfo& parent,
                                  int copyNo,
                                  bool isVisible,
                                  G4Colour color,
                                  bool forceSolid
                                  ){
    
    VolumeInfo info(name,offset,parent.centerInWorld);
    
    info.solid   = new G4Box( name, halfDim[0], halfDim[1], halfDim[2] );
    
    info.logical = new G4LogicalVolume( info.solid, material, name); 
    
    info.physical =  new G4PVPlacement( rot, 
                                        offset, 
                                        info.logical, 
                                        name, 
                                        parent.logical, 
                                        0, 
                                        copyNo,
                                        _config->getBool("g4.doSurfaceCheck",false));

    if ( isVisible ){

      // We need to manage the lifetime of the G4VisAttributes object.
      _visAttributes.push_back(G4VisAttributes(true, color));
      G4VisAttributes& visAtt = _visAttributes.back();

      // Finish the setting of visualization properties.
      visAtt.SetForceSolid(forceSolid);
      info.logical->SetVisAttributes(&visAtt);
    } 
    else{

      info.logical->SetVisAttributes(G4VisAttributes::Invisible);
    }

    // Save the volume information in case someone else needs to access it by name.
    addVolInfo(info);

    return info;
  }  // end of Mu2eWorld::nestBox

  // Create and place a G4Tubs inside a logical volume.
  VolumeInfo Mu2eWorld::nestTubs2 ( string const& name,
                                    double param[5],
                                    G4Material* material,
                                    G4RotationMatrix* rot,
                                    G4ThreeVector const& offset,
                                    const VolumeInfo& parent,
                                    int copyNo,
                                    bool isVisible,
                                    G4Colour color,
                                    bool forceSolid
                                    ){
    
    VolumeInfo info(name,offset,parent.centerInWorld);
    
    info.solid   = new G4Tubs( name, param[0], param[1], param[2], param[3], param[4]  );
    
    info.logical = new G4LogicalVolume( info.solid, material, name); 
    
    info.physical =  new G4PVPlacement( rot,
                                        offset,
                                        info.logical,
                                        name,
                                        parent.logical,
                                        0,
                                        copyNo,
                                        _config->getBool("g4.doSurfaceCheck",false));

    if ( isVisible ){

      // We need to manage the lifetime of the G4VisAttributes object.
      _visAttributes.push_back(G4VisAttributes(true, color));
      G4VisAttributes& visAtt = _visAttributes.back();

      // If I do not do this, then the rendering depends on what happens in
      // other parts of the code;  is there a G4 bug that causes something to be
      // unitialized?
      visAtt.SetForceAuxEdgeVisible(_config->getBool("g4.forceAuxEdgeVisible",false));

      // Finish the setting of visualization properties.
      visAtt.SetForceSolid(forceSolid);
      info.logical->SetVisAttributes(&visAtt);
    } 
    else{

      info.logical->SetVisAttributes(G4VisAttributes::Invisible);
    }

    // Save the volume information in case someone else needs to access it by name.
    addVolInfo(info);

    return info;
  }  // end of Mu2eWorld::nestTubs

  // Create and place a G4Cons inside a logical volume.
  VolumeInfo Mu2eWorld::nestCons2 ( string const& name,
                                    double param[7],
                                    G4Material* material,
                                    G4RotationMatrix* rot,
                                    G4ThreeVector const& offset,
                                    const VolumeInfo& parent,
                                    int copyNo,
                                    bool isVisible,
                                    G4Colour color,
                                    bool forceSolid
                                    ){
    
    VolumeInfo info(name,offset,parent.centerInWorld);
    
    info.solid   = new G4Cons( name, param[0], param[1], param[2], param[3], 
			             param[4], param[5], param[6]  );
    
    info.logical = new G4LogicalVolume( info.solid, material, name); 
    
    info.physical =  new G4PVPlacement( rot, offset, info.logical, name, parent.logical, 0, copyNo);

    if ( isVisible ){

      // We need to manage the lifetime of the G4VisAttributes object.
      _visAttributes.push_back(G4VisAttributes(true, color));
      G4VisAttributes& visAtt = _visAttributes.back();

      // If I do not do this, then the rendering depends on what happens in
      // other parts of the code;  is there a G4 bug that causes something to be
      // unitialized?
      visAtt.SetForceAuxEdgeVisible(_config->getBool("g4.forceAuxEdgeVisible",false));

      // Finish the setting of visualization properties.
      visAtt.SetForceSolid(forceSolid);
      info.logical->SetVisAttributes(&visAtt);
    } 
    else{

      info.logical->SetVisAttributes(G4VisAttributes::Invisible);
    }

    // Save the volume information in case someone else needs to access it by name.
    addVolInfo(info);

    return info;
  }  // end of Mu2eWorld::nestCons

  // Create and place a G4ExtrudedSolid inside a logical volume.
  VolumeInfo Mu2eWorld::nestExtrudedSolid2 ( string const& name,
					     double hz,
					     vector<double> & x,
					     vector<double> & y,
					     G4Material* material,
					     G4RotationMatrix* rot,
					     G4ThreeVector const& offset,
					     const VolumeInfo& parent,
					     int copyNo,
					     bool isVisible,
					     G4Colour color,
					     bool forceSolid
					     ){
    
    VolumeInfo info(name,offset,parent.centerInWorld);

    if( x.size()!=y.size() || x.size()==0 ) return info;

    vector<G4TwoVector> polygon;
    for( int i=0; i<x.size(); ++i ) polygon.push_back(G4TwoVector(x[i],y[i]));
    
    info.solid   = new G4ExtrudedSolid( name, polygon, hz,
					G4TwoVector(0.0,0.0), 1.0,
					G4TwoVector(0.0,0.0), 1.0 );
    
    info.logical = new G4LogicalVolume( info.solid, material, name); 
    
    info.physical =  new G4PVPlacement( rot, offset, info.logical, name, parent.logical, 0, copyNo);

    if ( isVisible ){

      // We need to manage the lifetime of the G4VisAttributes object.
      _visAttributes.push_back(G4VisAttributes(true, color));
      G4VisAttributes& visAtt = _visAttributes.back();

      // If I do not do this, then the rendering depends on what happens in
      // other parts of the code;  is there a G4 bug that causes something to be
      // unitialized?
      visAtt.SetForceAuxEdgeVisible(_config->getBool("g4.forceAuxEdgeVisible",false));

      // Finish the setting of visualization properties.
      visAtt.SetForceSolid(forceSolid);
      info.logical->SetVisAttributes(&visAtt);
    } 
    else{

      info.logical->SetVisAttributes(G4VisAttributes::Invisible);
    }

    // Save the volume information in case someone else needs to access it by name.
    addVolInfo(info);

    return info;
  }  // end of Mu2eWorld::nestCons


  //
  // Create and place a G4Torus inside a logical volume.
  // 
  VolumeInfo Mu2eWorld::nestTorus2 ( string const& name,
                                     double param[5],
                                     G4Material* material,
                                     G4RotationMatrix* rot,
                                     G4ThreeVector const& offset,
                                     const VolumeInfo& parent,
                                     int copyNo,
                                     bool isVisible,
                                     G4Colour color,
                                     bool forceSolid
                                     ){
    
    VolumeInfo info(name,offset,parent.centerInWorld);
    
    info.solid   = new G4Torus( name, param[0], param[1], param[2], param[3], param[4]  );
    
    info.logical = new G4LogicalVolume( info.solid, material, name); 
    
    info.physical =  new G4PVPlacement( rot,
                                        offset,
                                        info.logical,
                                        name,
                                        parent.logical,
                                        0,
                                        copyNo,
                                        _config->getBool("g4.doSurfaceCheck",false));
    
    if ( isVisible ){

      // We need to manage the lifetime of the G4VisAttributes object.
      _visAttributes.push_back(G4VisAttributes(true, color));
      G4VisAttributes& visAtt = _visAttributes.back();

      // If I do not do this, then the rendering depends on what happens in
      // other parts of the code;  is there a G4 bug that causes something to be
      // unitialized?
      visAtt.SetForceAuxEdgeVisible(_config->getBool("g4.forceAuxEdgeVisible",false));

      // Finish the setting of visualization properties.
      visAtt.SetForceSolid(forceSolid);
      info.logical->SetVisAttributes(&visAtt);
    } 
    else{

      info.logical->SetVisAttributes(G4VisAttributes::Invisible);
    }

    // Save the volume information in case someone else needs to access it by name.
    addVolInfo(info);

    return info;
  } // end of Mu2eWorld::nestTorus

  // Return the volume info mapped to the given key, throw if the key does not exist.
  // The syntax:
  //    VolumeInfo& xx =_volumeInfoList[key];
  // does not throw if the key does not exist.
  VolumeInfo& Mu2eWorld::locateVolInfo( const std::string key){
    std::map<std::string,VolumeInfo>::iterator i = _volumeInfoList.find(key);
    if ( i == _volumeInfoList.end() ){
      throw cms::Exception("GEOM")
        << "locateVolInfo cannot find the volume named: "
        << key 
        << "\n";
    }
    return i->second;
  } // end of Mu2eWorld::locateVolInfo

  // If the key already exists, throw. Otherwise add the (key, value) pair
  // to the map.
  // The syntax:
  //    VolumeInfo& xx;
  //   _volumeInfoList[key] = xx;
  // does not throw if the key already exists.
  void Mu2eWorld::addVolInfo( const VolumeInfo& info ){
    std::map<std::string,VolumeInfo>::iterator i = _volumeInfoList.find(info.name);
    if ( i != _volumeInfoList.end() ){
      throw cms::Exception("GEOM")
        << "locateVolInfo already has the key: "
        << info.name
        << "\n";
    }
    _volumeInfoList[info.name] = info;
  } // end of Mu2eWorld::addVolInfo


} // end namespace mu2e

