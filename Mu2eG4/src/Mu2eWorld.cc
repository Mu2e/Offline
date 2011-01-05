//
// Construct the Mu2e G4 world and serve information about that world.
//
// $Id: Mu2eWorld.cc,v 1.75 2011/01/05 21:04:47 genser Exp $
// $Author: genser $ 
// $Date: 2011/01/05 21:04:47 $
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

// C++ includes
#include <iostream>
#include <vector>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/constructWorldVolume.hh"
#include "Mu2eG4/inc/constructDirt.hh"
#include "Mu2eG4/inc/constructHall.hh"
#include "Mu2eG4/inc/constructProtonAbsorber.hh"
#include "Mu2eG4/inc/constructSteel.hh"
#include "Mu2eG4/inc/constructVirtualDetectors.hh"
#include "Mu2eG4/inc/constructDS.hh"
#include "Mu2eG4/inc/constructTS.hh"
#include "Mu2eG4/inc/constructPS.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/VirtualDetectorSD.hh"
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/CaloReadoutSD.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestTorus.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestCons.hh"
#include "Mu2eG4/inc/finishNesting.hh"
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
#include "TTrackerGeom/inc/TTracker.hh"

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
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
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


  // This is the callback called by G4 via G4VPhysicalVolume* WorldMaker::Construct()
  WorldInfo const* Mu2eWorld::construct(){

    _helper = &(*(edm::Service<G4Helper>()));

    // Get access to the master geometry system and its run time config.
    edm::Service<GeometryService> geom;
    _config = &(geom->config());

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
    CaloReadoutSD::setMu2eOriginInWorld( _mu2eOrigin );

    instantiateSensitiveDetectors();

    VolumeInfo worldVInfo = constructWorldVolume(_config);

    _info.worldPhys  = worldVInfo.physical;

    VolumeInfo dirtInfo  = constructDirt( worldVInfo,_config );
    VolumeInfo hallInfo  = constructHall( dirtInfo,_config );

    // Define the hall origin in Mu2e coordinates.
    _hallOriginInMu2e = hallInfo.centerInWorld - _mu2eOrigin;

    constructDS(hallInfo,_config);
    constructTS(hallInfo,_config);
    constructPS(hallInfo,_config, _primaryProtonGunOrigin, _primaryProtonGunRotation);// input/output params

    VolumeInfo trackerInfo = constructTracker();
    VolumeInfo targetInfo  = constructTarget();

    constructProtonAbsorber(_config);

    constructVirtualDetectors(_config);

    // These are just placeholders for now - and might be misnamed.
    constructCal();
    constructMagnetYoke();
    constructCRV();

    constructSteel(hallInfo,_config);

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


  // Choose the selected tracker and build it.
  VolumeInfo Mu2eWorld::constructTracker(){

    // The tracker is built inside this volume.
    VolumeInfo const & detSolDownstreamVacInfo = _helper->locateVolInfo("ToyDS3Vacuum");

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
      // Hack alert: These belong in constructTracker
      trackerInfo.name = "TrackerMother"; // this belongs to construct..., some of them do it now
      _helper->addVolInfo(trackerInfo);
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
    VolumeInfo const & detSolUpstreamVacInfo   = _helper->locateVolInfo("ToyDS2Vacuum");

    // z Position of the center of the DS solenoid parts, given in the Mu2e coordinate system.
    double z0DSup   = detSolUpstreamVacInfo.centerInWorld.z()+_hallOriginInMu2e.z();

    // Buid the stopping target
    VolumeInfo targetInfo = ( _config->getBool("hasTarget",false) ) ? 

      constructStoppingTarget( detSolUpstreamVacInfo.logical, 
                                            z0DSup,
                                            *_config )
      :

      constructDummyStoppingTarget( detSolUpstreamVacInfo.logical, 
                                                 z0DSup,
                                                 *_config );
    return targetInfo;

  } // end Mu2eWorld::constructTarget

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
    G4LogicalVolume* ds2Vacuum = _helper->locateVolInfo("ToyDS2Vacuum").logical;
    G4LogicalVolume* ds3Vacuum = _helper->locateVolInfo("ToyDS3Vacuum").logical;

    vector<G4LogicalVolume*> psVacua;
    psVacua.push_back( _helper->locateVolInfo("PS1Vacuum").logical );

    vector<G4LogicalVolume*> tsVacua;
    tsVacua.push_back( _helper->locateVolInfo("ToyTS1Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("ToyTS2Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("ToyTS3Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("ToyTS4Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("ToyTS5Vacuum").logical );

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
    G4LogicalVolume* ds2Vacuum = _helper->locateVolInfo("ToyDS2Vacuum").logical;
    G4LogicalVolume* ds3Vacuum = _helper->locateVolInfo("ToyDS3Vacuum").logical;

    vector<G4LogicalVolume*> psVacua;
    psVacua.push_back( _helper->locateVolInfo("PS1Vacuum").logical );

    vector<G4LogicalVolume*> tsVacua;
    tsVacua.push_back( _helper->locateVolInfo("ToyTS1Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("ToyTS2Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("ToyTS3Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("ToyTS4Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("ToyTS5Vacuum").logical );

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

    G4LogicalVolume* ds2Vacuum      = _helper->locateVolInfo("ToyDS2Vacuum").logical;
    G4LogicalVolume* ds3Vacuum      = _helper->locateVolInfo("ToyDS3Vacuum").logical;
    G4LogicalVolume* tracker        = _helper->locateVolInfo("TrackerMother").logical;
    G4LogicalVolume* stoppingtarget = _helper->locateVolInfo("StoppingTargetMother").logical;

    vector<G4LogicalVolume*> psVacua;
    psVacua.push_back( _helper->locateVolInfo("PS1Vacuum").logical );

    vector<G4LogicalVolume*> tsVacua;
    tsVacua.push_back( _helper->locateVolInfo("ToyTS1Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("ToyTS2Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("ToyTS3Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("ToyTS4Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("ToyTS5Vacuum").logical );

    // We may make separate G4UserLimits objects per logical volume but we choose not to.
    //_stepLimits.push_back( G4UserLimits(maxStep) );
    //G4UserLimits* stepLimit = &(_stepLimits.back());

    AntiLeakRegistry& reg = edm::Service<G4Helper>()->antiLeakRegistry();
    G4UserLimits* stepLimit = reg.add( G4UserLimits(maxStep) );
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

    VolumeInfo const & detSolDownstreamVacInfo = _helper->locateVolInfo("ToyDS3Vacuum");

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

  // instantiateSensitiveDetectors

  void Mu2eWorld::instantiateSensitiveDetectors(){

    G4SDManager* SDman      = G4SDManager::GetSDMpointer();

    // G4 takes ownership and will delete the detectors at the job end

    StrawSD* strawSD        = new StrawSD(          SensitiveDetectorName::StrawGasVolume(),  *_config);
    //strawSD->SetVerboseLevel(1);
    SDman->AddNewDetector(strawSD); 

    VirtualDetectorSD* vdSD = new VirtualDetectorSD(SensitiveDetectorName::VirtualDetector(), *_config);
    SDman->AddNewDetector(vdSD);

    CaloCrystalSD* ccSD     = new CaloCrystalSD(    SensitiveDetectorName::CaloCrystal(),     *_config);
    SDman->AddNewDetector(ccSD);

    CaloReadoutSD* crSD     = new CaloReadoutSD(    SensitiveDetectorName::CaloReadout(),     *_config);
    SDman->AddNewDetector(crSD);

  } // instantiateSensitiveDetectors

} // end namespace mu2e
