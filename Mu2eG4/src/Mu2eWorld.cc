//
// Construct the Mu2e G4 world and serve information about that world.
//
// $Id: Mu2eWorld.cc,v 1.92 2011/05/26 22:11:19 genser Exp $
// $Author: genser $
// $Date: 2011/05/26 22:11:19 $
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
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/constructWorldVolume.hh"
#include "Mu2eG4/inc/constructDirt.hh"
#include "Mu2eG4/inc/constructHall.hh"
#include "Mu2eG4/inc/constructProtonAbsorber.hh"
#include "Mu2eG4/inc/constructSteel.hh"
#include "Mu2eG4/inc/constructCRV.hh"
#include "Mu2eG4/inc/constructNeutronAbsorber.hh"
#include "Mu2eG4/inc/constructMBS.hh"
#include "Mu2eG4/inc/constructVirtualDetectors.hh"
#include "Mu2eG4/inc/constructDS.hh"
#include "Mu2eG4/inc/constructTS.hh"
#include "Mu2eG4/inc/constructPS.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/StrawSD.hh"
//#include "Mu2eG4/inc/ITGasLayerSD.hh"
#include "Mu2eG4/inc/ITGasLayerSD_Hexagonal.hh"
#include "Mu2eG4/inc/ITGasLayerSD_Square.hh"
#include "Mu2eG4/inc/VirtualDetectorSD.hh"
#include "Mu2eG4/inc/StoppingTargetSD.hh"
#include "Mu2eG4/inc/CRSScintillatorBarSD.hh"
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

    _helper = &(*(art::ServiceHandle<G4Helper>()));

    // Get access to the master geometry system and its run time config.
    art::ServiceHandle<GeometryService> geom;
    _config = &(geom->config());

    // Construct all of the Mu2e world, hall, detectors, beamline ...
    constructWorld();

    return &_info;
  }

  // Construct all of the Mu2e world, hall, detectors, beamline ...
  void Mu2eWorld::constructWorld(){

    int static const diagLevel = 0;

    // If you play with the order of these calls, you may break things.
    defineMu2eOrigin();
    VolumeInfo::setMu2eOriginInWorld( _mu2eOrigin );
    VirtualDetectorSD::setMu2eOriginInWorld( _mu2eOrigin );
    StoppingTargetSD::setMu2eOriginInWorld( _mu2eOrigin );
    CaloCrystalSD::setMu2eOriginInWorld( _mu2eOrigin );
    CaloReadoutSD::setMu2eOriginInWorld( _mu2eOrigin );
    CRSScintillatorBarSD::setMu2eOriginInWorld( _mu2eOrigin );
    if ( _config->getBool("hasITracker",false) ) {
            ITGasLayerSD::setMu2eDetCenterInWorld( _mu2eDetectorOrigin -
                            G4ThreeVector(0.0,0.0,12000-_config->getDouble("itracker.z0",0.0)) );
    }

    instantiateSensitiveDetectors();

    VolumeInfo worldVInfo = constructWorldVolume(_config);

    _info.worldPhys  = worldVInfo.physical;

    if ( diagLevel > 0) {
      cout << __func__ << " worldVInfo.centerInParent : " <<  worldVInfo.centerInParent << endl;
      cout << __func__ << " worldVInfo.centerInWorld  : " <<  worldVInfo.centerInWorld  << endl;
    }

    VolumeInfo dirtInfo  = constructDirt( worldVInfo,_config );

    if ( diagLevel > 0) {
      cout << __func__ << " dirtInfo.centerInParent   : " <<  dirtInfo.centerInParent << endl;
      cout << __func__ << " dirtInfo.centerInWorld    : " <<  dirtInfo.centerInWorld  << endl;
    }

    VolumeInfo hallInfo  = constructHall( dirtInfo,_config );

    if ( diagLevel > 0) {
      cout << __func__ << " hallInfo.centerInParent   : " <<  hallInfo.centerInParent << endl;
      cout << __func__ << " hallInfo.centerInWorld    : " <<  hallInfo.centerInWorld  << endl;
    }

    // Define the hall origin in Mu2e coordinates.
    _hallOriginInMu2e = hallInfo.centerInWorld - _mu2eOrigin;

    constructDS(hallInfo,_config);
    constructTS(hallInfo,_config);
    constructPS(hallInfo,_config, _primaryProtonGunOrigin, _primaryProtonGunRotation);// input/output params

    VolumeInfo trackerInfo = constructTracker();
    VolumeInfo targetInfo  = constructTarget();

    constructProtonAbsorber(_config);

    // These are just placeholders for now - and might be misnamed.
    constructCal();
    constructMagnetYoke();

    if ( _config->getBool("hasCosmicRayShield",false) ) {
      constructSteel(hallInfo,_config);
      constructCRV(hallInfo,_config);
    }

    if ( _config->getBool("hasNeutronAbsorber",false) ) {
      constructNeutronAbsorber(_config);
    }

    if ( _config->getBool("hasMBS",false) ) {
      constructMBS(_config);
    }

    constructVirtualDetectors(_config); // beware of the placement order of this function

    mf::LogInfo log("GEOM");
    log << "Mu2e Origin:          " << _mu2eOrigin           << "\n";
    log << "Mu2e Detector Origin: " << _mu2eDetectorOrigin   << "\n";
    log << "Cosmic Ref:           " << _cosmicReferencePoint << "\n";
    log << "Hall Origin in Mu2e:  " << _hallOriginInMu2e     << "\n";

    // Create magnetic fields and managers only after all volumes have been defined.
    constructBFieldAndManagers();
    constructStepLimiters();
    if ( _config->getBool("hasITracker",false) ) {
            constructITStepLimiters();
    }

  }

  // Convert to base units for all of the items in the vector.
  void Mu2eWorld::setUnits( vector<double>& V, G4double unit ){
    for ( vector<double>::iterator b=V.begin(), e=V.end();
          b!=e; ++b){
      *b *= unit;
    }
  }

  void Mu2eWorld::defineMu2eOrigin(){

    int static const diagLevel = 0;

    // Dimensions of the world.
    vector<double> worldHLen;
    _config->getVectorDouble("world.halfLengths", worldHLen, 3);

    // Floor thickness.
    double floorThick = _config->getDouble("hall.floorThick");

    // Top of the floor in G4 world coordinates.
    double yFloor = -worldHLen[1] + floorThick;

    if ( diagLevel > 0) {
      cout << __func__ << " yFlor : " <<  yFloor  << endl;
    }

    // The height above the floor of the y origin of the Mu2e coordinate system.
    double yOriginHeight = _config->getDouble("world.mu2eOrigin.height" )*CLHEP::mm;

    // Position of the origin of the mu2e coordinate system
    _mu2eOrigin = G4ThreeVector(
                                _config->getDouble("world.mu2eOrigin.xoffset")*CLHEP::mm,
                                yFloor + yOriginHeight,
                                _config->getDouble("world.mu2eOrigin.zoffset")*CLHEP::mm
                                );

    if ( diagLevel > 0) {
      cout << __func__ << " _mu2eOrigin : " <<  _mu2eOrigin  << endl;
    }

    // Origin used to construct the MECO detector.
    // Magic number to fix:
    _mu2eDetectorOrigin = _mu2eOrigin + G4ThreeVector( -3904., 0., 12000.);

    if ( diagLevel > 0) {
      cout << __func__ << " _mu2eDetectorOrigin : " <<  _mu2eDetectorOrigin  << endl;
    }

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

    _dirtG4Ymax = ySurface;
    _dirtG4Ymin = yCeilingOutside;

    if ( diagLevel > 0) {
      cout << __func__ << " yEverest : " <<  yEverest  << endl;
    }

    // Selfconsistency check.
    if ( yEverest > 2.*worldHLen[1] ){
      throw cet::exception("GEOM")
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
    // is it better to use here detSolUpstreamVacInfo.centerInMu2e().z() ?
    double z0DSup   = detSolUpstreamVacInfo.centerInWorld.z()+_hallOriginInMu2e.z();

    cout << "_mu2eOrigin.z()=" << _mu2eOrigin.z() << endl;
    cout << "detSolUpstreamVacInfo.centerInWorld.z()=" << detSolUpstreamVacInfo.centerInWorld.z() << endl;
    cout << "detSolUpstreamVacInfo.centerInMu2e.z()=" << detSolUpstreamVacInfo.centerInMu2e().z() << endl;
    cout << "_hallOriginInMu2e.z()=" << _hallOriginInMu2e.z() << endl;

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
  // the relevant volumes or to the world

  void Mu2eWorld::constructBFieldAndManagers(){

    // Figure out which magnetic field managers are needed.
    int dsFieldForm    = _config->getInt("detSolFieldForm", dsModelUniform);

    // Decide on the G4 Stepper

    bool needDSUniform = (dsFieldForm == dsModelSplit || dsFieldForm == dsModelUniform );

    string stepper = _config->getString("g4.stepper","G4SimpleRunge");

    // Create field manager for the uniform DS field.
    if ( needDSUniform){
      // Handle to the BField manager.
      GeomHandle<BFieldManager> bfMgr;
      _dsUniform = FieldMgr::forUniformField( bfMgr->getDSUniformValue()*CLHEP::tesla, _mu2eOrigin );
    }

    // Create global field managers; don't use FieldMgr here to avoid problem with ownership

    G4MagneticField * _field = new DSField("", _mu2eOrigin);
    G4Mag_UsualEqRhs * _rhs  = new G4Mag_UsualEqRhs(_field);
    G4MagIntegratorStepper * _stepper;
    if ( stepper  == "G4ClassicalRK4" ) {
      _stepper = new G4ClassicalRK4(_rhs);
    } else if ( stepper  == "G4ImplicitEuler" ) {
      _stepper = new G4ImplicitEuler(_rhs);
    } else if ( stepper  == "G4ExplicitEuler" ) {
      _stepper = new G4ExplicitEuler(_rhs);
    } else if ( stepper  == "G4SimpleHeum" ) {
      _stepper = new G4SimpleHeum(_rhs);
    } else if ( stepper  == "G4HelixImplicitEuler" ) {
      _stepper = new G4HelixImplicitEuler(_rhs);
    } else if ( stepper  == "G4HelixSimpleRunge" ) {
      _stepper = new G4HelixSimpleRunge(_rhs);
    } else {
      _stepper = new G4SimpleRunge(_rhs);
    }
    G4ChordFinder * _chordFinder = new G4ChordFinder(_field,1.0e-2*CLHEP::mm,_stepper);
    G4FieldManager * _manager = new G4FieldManager(_field,_chordFinder,true);

    // G4TransportationManager takes ownership of _manager
    G4TransportationManager::GetTransportationManager()->SetFieldManager(_manager);

    // Define uniform field region in the detector solenoid, if neccessary

    G4LogicalVolume* ds2Vacuum = _helper->locateVolInfo("ToyDS2Vacuum").logical;
    G4LogicalVolume* ds3Vacuum = _helper->locateVolInfo("ToyDS3Vacuum").logical;

    if (dsFieldForm == dsModelUniform  ){
      ds2Vacuum->SetFieldManager( _dsUniform->manager(), true);
      ds3Vacuum->SetFieldManager( _dsUniform->manager(), true);
    } else if ( dsFieldForm == dsModelSplit ){
      ds3Vacuum->SetFieldManager( _dsUniform->manager(), true);
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

    _manager->SetDeltaOneStep(deltaOneStep);
    _manager->SetDeltaIntersection(newUpstreamDeltaI);
    _chordFinder->SetDeltaChord(deltaChord);

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

    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
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

  void Mu2eWorld::constructITStepLimiters(){

    bool physicalStep =  _config->getBool("itracker.usePhysicalStep",false);
    // Maximum step length, in mm.
    double maxStep = 10.0;
    if (physicalStep){
            maxStep = 10.0/12.0;
    }else {
            maxStep = _config->getDouble("itracker.freePath", 0.5);
    }
    G4LogicalVolume* tracker        = _helper->locateVolInfo("TrackerMother").logical;

    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
    G4UserLimits* stepLimit = reg.add( G4UserLimits(maxStep) );

    GeomHandle<ITracker> itracker;
    SuperLayer *iSpl;
    int ilay;
    for (int ispls=0; ispls<itracker->nSuperLayers(); ispls++){
            iSpl = itracker->getSuperLayer(ispls);
            for ( ilay=0; ilay<iSpl->nLayers(); ilay++){
                    iSpl->getLayer(ilay)->getDetail();
            }

    }
    G4VPhysicalVolume *iDau;
    //cout<<"N IT daughter: "<<tracker->TotalVolumeEntities()<<endl;
    for (int iDaughter=0; iDaughter<tracker->TotalVolumeEntities(); iDaughter++){
            iDau = tracker->GetDaughter(iDaughter);
            if (!iDau) break;
            //cout<<"Vol Name "<< iDau->GetName()<<" is Tracking: "<<iDau->GetName().contains("volS")<<endl;
            if (iDau->GetName().contains("volS")) iDau->GetLogicalVolume()->SetUserLimits(stepLimit);
    }

    cout<<"IT Step limits set"<<endl;


  } // end Mu2eWorld::constructITStepLimiters(){


  // Construct calorimeter if needed.
  void Mu2eWorld::constructCal(){

    if ( ! _config->getBool("hasCalorimeter",false) ) return;

    VolumeInfo const & detSolDownstreamVacInfo = _helper->locateVolInfo("ToyDS3Vacuum");

    double z0DSdown = detSolDownstreamVacInfo.centerInWorld.z()+_hallOriginInMu2e.z();

    constructCalorimeter( detSolDownstreamVacInfo,
                          -z0DSdown,
                          *_config );
  }


  // A place holder for now.
  void Mu2eWorld::constructMagnetYoke(){
  } // end Mu2eWorld::constructMagnetYoke

  // instantiateSensitiveDetectors

  void Mu2eWorld::instantiateSensitiveDetectors(){

    G4SDManager* SDman      = G4SDManager::GetSDMpointer();

    // G4 takes ownership and will delete the detectors at the job end

    if ( _config->getBool("hasITracker",false) ) {
      GeomHandle<ITracker> itracker;
      ITGasLayerSD* itrackerSD=0x0;
      if ( itracker->geomType()==ITracker::Hexagonal )
        itrackerSD = new ITGasLayerSD_Hexagonal(SensitiveDetectorName::ItrackerGasVolume(), *_config);
      else if ( itracker->geomType()==ITracker::Square )
        itrackerSD = new ITGasLayerSD_Square(   SensitiveDetectorName::ItrackerGasVolume(), *_config);
      //itrackerSD->SetVerboseLevel(1);
      SDman->AddNewDetector(itrackerSD);
    }
    else {
      StrawSD* strawSD      = new StrawSD(  SensitiveDetectorName::StrawGasVolume(),  *_config);
      //strawSD->SetVerboseLevel(1);
      SDman->AddNewDetector(strawSD);
    }

    VirtualDetectorSD* vdSD = new VirtualDetectorSD(SensitiveDetectorName::VirtualDetector(), *_config);
    SDman->AddNewDetector(vdSD);

    CaloCrystalSD* ccSD     = new CaloCrystalSD(    SensitiveDetectorName::CaloCrystal(),     *_config);
    SDman->AddNewDetector(ccSD);

    CaloReadoutSD* crSD     = new CaloReadoutSD(    SensitiveDetectorName::CaloReadout(),     *_config);
    SDman->AddNewDetector(crSD);

    StoppingTargetSD* stSD = new StoppingTargetSD(  SensitiveDetectorName::StoppingTarget(),  *_config);
    SDman->AddNewDetector(stSD);

    CRSScintillatorBarSD* sbSD = new CRSScintillatorBarSD(SensitiveDetectorName::CRSScintillatorBar(), *_config);
    SDman->AddNewDetector(sbSD);


  } // instantiateSensitiveDetectors

} // end namespace mu2e
