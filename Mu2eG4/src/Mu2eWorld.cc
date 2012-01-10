//
// Construct the Mu2e G4 world and serve information about that world.
//
// $Id: Mu2eWorld.cc,v 1.112 2012/01/10 22:28:43 mu2ecvs Exp $
// $Author: mu2ecvs $
// $Date: 2012/01/10 22:28:43 $
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
#include "Mu2eG4/inc/constructProtonBeamDump.hh"
#include "Mu2eG4/inc/constructProtonAbsorber.hh"
#include "Mu2eG4/inc/constructSteel.hh"
#include "Mu2eG4/inc/constructCRV.hh"
#include "Mu2eG4/inc/constructExtMonUCI.hh"
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
#include "Mu2eG4/inc/ExtMonFNAL_SD.hh"
#include "Mu2eG4/inc/ExtMonUCITofSD.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestTorus.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestCons.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/ITrackerBuilder.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/Mu2eBuilding.hh"
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
#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCI.hh"

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
#include "G4GDMLParser.hh"

#include "Mu2eG4/inc/DSField.hh"
#include "Mu2eG4/inc/FieldMgr.hh"

using namespace std;

namespace mu2e {

  Mu2eWorld::Mu2eWorld():
    _helper(0)
  {
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
  G4VPhysicalVolume * Mu2eWorld::construct(){
    
    _helper = &(*(art::ServiceHandle<G4Helper>()));

    // Get access to the master geometry system and its run time config.
    art::ServiceHandle<GeometryService> geom;
    _config = &(geom->config());

    // Construct all of the Mu2e world, hall, detectors, beamline ...
    return constructWorld();
  }
  
  // Construct all of the Mu2e world, hall, detectors, beamline ...
  G4VPhysicalVolume * Mu2eWorld::constructWorld(){
    
    int static const diagLevel = _config->getInt("world.verbosityLevel", 0);
    art::ServiceHandle<GeometryService> geom;

    // If you play with the order of these calls, you may break things.
    GeomHandle<WorldG4> worldGeom;
    if ( _config->getBool("hasITracker",false) ) {
      ITGasLayerSD::setMu2eDetCenterInWorld( 
					    GeomHandle<Mu2eBuilding>()->trackerOriginInMu2e() + worldGeom->mu2eOriginInWorld() 
					    - G4ThreeVector(0.0,0.0,12000-_config->getDouble("itracker.z0",0.0)) );
    }

    instantiateSensitiveDetectors();

    VolumeInfo worldVInfo = constructWorldVolume(_config);

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
      cout << __func__ << " hallInfo.centerInMu2e()   : " <<  hallInfo.centerInMu2e() << endl;
    }

    constructProtonBeamDump(hallInfo, *_config);

    constructDS(hallInfo,_config);
    constructTS(hallInfo,_config);
    constructPS(hallInfo,_config);

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

    if ( _config->getBool("hasExtMonUCI",false) ) {
      constructExtMonUCI(hallInfo, *_config);
    }

    constructVirtualDetectors(_config); // beware of the placement order of this function

    mf::LogInfo log("GEOM");
    log << "Mu2e Origin:          " << worldGeom->mu2eOriginInWorld() << "\n";
    log << "Mu2e Detector Origin: " << GeomHandle<Mu2eBuilding>()->trackerOriginInMu2e() + worldGeom->mu2eOriginInWorld()   << "\n";
    log << "Cosmic Ref:           " << worldGeom->cosmicReferencePoint() << "\n";

    // Create magnetic fields and managers only after all volumes have been defined.
    constructBFieldAndManagers();
    constructStepLimiters();
    if ( _config->getBool("hasITracker",false) ) {
            constructITStepLimiters();
    }

    // Write out mu2e geometry into a gdml file.
    if ( _config->getBool("writeGDML",false) ) {
      string gdmlFileName = _config->getString("GDMLFileName","mu2e.gdml");
      G4GDMLParser parser;
      parser.Write(gdmlFileName, worldVInfo.logical);
    }

    return worldVInfo.physical;
  }

  // Convert to base units for all of the items in the vector.
  void Mu2eWorld::setUnits( vector<double>& V, G4double unit ){
    for ( vector<double>::iterator b=V.begin(), e=V.end();
          b!=e; ++b){
      *b *= unit;
    }
  }

  // Choose the selected tracker and build it.
  VolumeInfo Mu2eWorld::constructTracker(){

    // The tracker is built inside this volume.
    VolumeInfo const & detSolDownstreamVacInfo = _helper->locateVolInfo("ToyDS3Vacuum");

    // z Position of the center of the DS solenoid parts, given in the Mu2e coordinate system.
    double z0DSdown = detSolDownstreamVacInfo.centerInMu2e().z();

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
        trackerInfo = constructTTrackerv3( detSolDownstreamVacInfo, z0DSdown, *_config );
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

    cout << "detSolUpstreamVacInfo.centerInWorld.z()=" << detSolUpstreamVacInfo.centerInWorld.z() << endl;
    cout << "detSolUpstreamVacInfo.centerInMu2e().z() =" << detSolUpstreamVacInfo.centerInMu2e().z() << endl;

    // Buid the stopping target
    VolumeInfo targetInfo = ( _config->getBool("hasTarget",false) ) ?

      constructStoppingTarget( detSolUpstreamVacInfo,
                               *_config )
      :

      constructDummyStoppingTarget( detSolUpstreamVacInfo,
                                    *_config );
    return targetInfo;

  } // end Mu2eWorld::constructTarget

  // Construct the magnetic field managers and attach them to
  // the relevant volumes or to the world

  void Mu2eWorld::constructBFieldAndManagers(){
    
    GeomHandle<WorldG4> worldGeom;

    // Get some information needed further

    VolumeInfo const & ds2VacuumVacInfo = _helper->locateVolInfo("ToyDS2Vacuum");
    VolumeInfo const & ds3VacuumVacInfo = _helper->locateVolInfo("ToyDS3Vacuum");
    G4LogicalVolume* ds2Vacuum = ds2VacuumVacInfo.logical;
    G4LogicalVolume* ds3Vacuum = ds3VacuumVacInfo.logical;

    double ds2HL = static_cast<G4Tubs*>(ds2VacuumVacInfo.solid)->GetZHalfLength();
    G4ThreeVector ds2Z0 = ds2VacuumVacInfo.centerInWorld;
    G4ThreeVector beamZ0( ds2Z0.x(), ds2Z0.y(), ds2Z0.z()+ds2HL );

    // Figure out which magnetic field managers are needed.
    int dsFieldForm    = _config->getInt("detSolFieldForm", dsModelUniform);

    // Decide on the G4 Stepper

    bool needDSUniform = (dsFieldForm == dsModelSplit || dsFieldForm == dsModelUniform );
    bool needDSGradient = false;

    string stepper = _config->getString("g4.stepper","G4SimpleRunge");

    // Create field manager for the uniform DS field.
    if ( needDSUniform){
      // Handle to the BField manager.
      GeomHandle<BFieldManager> bfMgr;
      _dsUniform = FieldMgr::forUniformField( bfMgr->getDSUniformValue()*CLHEP::tesla, worldGeom->mu2eOriginInWorld() );
      
      // Create field manager for the gradient field in DS3
      if( fabs(bfMgr->getDSGradientValue().z())>1.0e-9 ) {
	needDSGradient = true;
	_dsGradient = FieldMgr::forGradientField( bfMgr->getDSUniformValue().z()*CLHEP::tesla,
						  bfMgr->getDSGradientValue().z()*CLHEP::tesla/CLHEP::m,
						  beamZ0 );
      }
    }

    // Create global field managers; don't use FieldMgr here to avoid problem with ownership

    G4MagneticField * _field = new DSField("", worldGeom->mu2eOriginInWorld());
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

    if (dsFieldForm == dsModelUniform  ){
      ds2Vacuum->SetFieldManager( _dsUniform->manager(), true);
      cout << "Use uniform field in DS2" << endl;
    }
    if ( dsFieldForm == dsModelUniform || dsFieldForm == dsModelSplit ){
      if( needDSGradient ) {
	ds3Vacuum->SetFieldManager( _dsGradient->manager(), true);
	cout << "Use gradient field in DS3" << endl;
      } else {
	ds3Vacuum->SetFieldManager( _dsUniform->manager(), true);
	cout << "Use uniform field in DS3" << endl;
      }
      /*
      if( _config->getBool("hasMBS",false) ) {
	VolumeInfo const & MBSMotherInfo = _helper->locateVolInfo("MBSMother");
	G4LogicalVolume* mbsMother = MBSMotherInfo.logical;
	mbsMother->SetFieldManager(0,true);
      }
      */
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

    if ( _dsGradient.get() != 0 ){
      G4double deltaIntersection = 0.00001*CLHEP::mm;
      _dsGradient->manager()->SetDeltaIntersection(deltaIntersection);
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

    double z0DSdown = detSolDownstreamVacInfo.centerInMu2e().z();

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

    ExtMonFNAL_SD* emfSD     = new ExtMonFNAL_SD(    SensitiveDetectorName::ExtMonFNAL(),     *_config);
    SDman->AddNewDetector(emfSD);

    ExtMonUCITofSD* emuTofSD     = new ExtMonUCITofSD(    SensitiveDetectorName::ExtMonUCITof(),     *_config);
    SDman->AddNewDetector(emuTofSD);

    StoppingTargetSD* stSD = new StoppingTargetSD(  SensitiveDetectorName::StoppingTarget(),  *_config);
    SDman->AddNewDetector(stSD);

    CRSScintillatorBarSD* sbSD = new CRSScintillatorBarSD(SensitiveDetectorName::CRSScintillatorBar(), *_config);
    SDman->AddNewDetector(sbSD);


  } // instantiateSensitiveDetectors

} // end namespace mu2e
