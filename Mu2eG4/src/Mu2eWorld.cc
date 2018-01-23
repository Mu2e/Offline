//
// Construct the Mu2e G4 world and serve information about that world.
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
#include <string>
#include <vector>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "ConfigTools/inc/checkForStale.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/constructWorldVolume.hh"
#include "Mu2eG4/inc/constructHall.hh"
#include "Mu2eG4/inc/constructProtonBeamDump.hh"
#include "Mu2eG4/inc/constructProtonAbsorber.hh"
#include "Mu2eG4/inc/constructCRV.hh"
#include "Mu2eG4/inc/constructExternalShielding.hh"
#include "Mu2eG4/inc/constructServicesGeom.hh"
#include "Mu2eG4/inc/constructTSdA.hh"
#include "Mu2eG4/inc/constructMBS.hh"
#include "Mu2eG4/inc/constructVirtualDetectors.hh"
#include "Mu2eG4/inc/constructVisualizationRegions.hh"
#include "Mu2eG4/inc/constructDS.hh"
#include "Mu2eG4/inc/constructMSTM.hh"
#include "Mu2eG4/inc/constructSTM.hh"
#include "Mu2eG4/inc/constructTS.hh"
#include "Mu2eG4/inc/constructPS.hh"
#include "Mu2eG4/inc/constructPSEnclosure.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/CaloReadoutSD.hh"
#include "Mu2eG4/inc/CaloReadoutCardSD.hh"
#include "Mu2eG4/inc/CaloCrateSD.hh"
#include "Mu2eG4/inc/ExtMonFNALPixelSD.hh"
#include "Mu2eG4/inc/TrackerWireSD.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/TTrackerPlaneSupportSD.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestTorus.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestCons.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/TransportSolenoid.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "Mu2eG4/inc/constructTTracker.hh"
#include "Mu2eG4/inc/constructDummyTracker.hh"
#include "Mu2eG4/inc/constructStoppingTarget.hh"
#include "Mu2eG4/inc/constructDummyStoppingTarget.hh"
#include "Mu2eG4/inc/constructDiskCalorimeter.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

// G4 includes
#include "G4Threading.hh"
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
#include "G4Mag_SpinEqRhs.hh"
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

#include "Mu2eG4/inc/Mu2eGlobalField.hh"

#include "boost/regex.hpp"

using namespace std;

namespace mu2e {

  Mu2eWorld::Mu2eWorld(SensitiveDetectorHelper *sdHelper/*no ownership passing*/)
    : sdHelper_(sdHelper)
    , pset_(fhicl::ParameterSet())

    , activeWr_Wl_SD_(_config.getBool("ttracker.ActiveWr_Wl_SD",false))
    , writeGDML_(_config.getBool("writeGDML",false))
    , gdmlFileName_(_config.getString("GDMLFileName","mu2e.gdml"))
    , g4stepperName_(_config.getString("g4.stepper","G4SimpleRunge"))
    , bfieldMaxStep_(_config.getDouble("bfield.maxStep", 20.))
    , worldVolume_(NULL)
  {
    _verbosityLevel = _config.getInt("world.verbosityLevel", 0);
  }

  Mu2eWorld::Mu2eWorld(const fhicl::ParameterSet& pset,
                       SensitiveDetectorHelper *sdHelper//,/*no ownership passing*/
                       //BFieldConfig *bfieldconf
                       )
    : sdHelper_(sdHelper)
    , pset_(pset)
//    , bfConfig(bfieldconf)

    , activeWr_Wl_SD_(true)
    , writeGDML_(pset.get<bool>("debug.writeGDML"))
    , gdmlFileName_(pset.get<std::string>("debug.GDMLFileName"))
    , g4stepperName_(pset.get<std::string>("physics.stepper"))
    , bfieldMaxStep_(pset.get<double>("physics.bfieldMaxStep"))
    , worldVolume_(NULL)
  {
    _verbosityLevel = pset.get<int>("debug.worldVerbosityLevel");
  }


    // This is the callback called by G4 via G4VPhysicalVolume* WorldMaker::Construct()
  G4VPhysicalVolume * Mu2eWorld::construct(){
        // Construct all of the Mu2e world, hall, detectors, beamline ...
      
      if (G4Threading::G4GetThreadId() <= 0){
          std::cout << "We are in Mu2eWorld::construct()" << std::endl;
      }

      
      
        return constructWorld();
  }
    
    // This is the callback called by G4 via void WorldMaker::ConstructSDandField()
  void Mu2eWorld::constructSDandField(){
      
      
      instantiateSensitiveDetectors();
      
      if (G4Threading::G4GetThreadId() == 0){
          std::cout << "We are in Mu2eWorld::constructSDandField() in Thread 0" << std::endl;
          G4SDManager* SDman = G4SDManager::GetSDMpointer();
          SDman->ListTree();
          
      }
      
      if (G4Threading::IsMasterThread() == true){
          std::cout << "We are in Mu2eWorld::constructSDandField() in the master thread" << std::endl;
          G4SDManager* SDman = G4SDManager::GetSDMpointer();
          SDman->ListTree();
          
      }

        //constructWorldSD();
      
      constructBFieldAndManagers();
      
  }

    //called by constructSDandField() to construct Sensitive Detectors
  void Mu2eWorld::constructWorldSD(){
  }

    
//*********************************************************************
//HERE'S the NEW FIELD FUNCTION
  void Mu2eWorld::constructWorldField(){

          //constructBFieldAndManagers();
    
  }
//END NEW FIELD FUNCTION
//*********************************************************************
    

  // Construct all of the Mu2e world, hall, detectors, beamline ...
  G4VPhysicalVolume * Mu2eWorld::constructWorld(){

      
      // If you play with the order of these calls, you may break things.
      GeomHandle<WorldG4> worldGeom;
      G4ThreeVector tmpTrackercenter = GeomHandle<DetectorSystem>()->getOrigin();

      if (activeWr_Wl_SD_) {
          TrackerWireSD::setMu2eDetCenterInWorld( tmpTrackercenter );
      }

//    instantiateSensitiveDetectors(); MOVED!!!!!!!!!!!!!
      
      
      // ************************  BEGIN: NO CHANGES NEEDED   ************************
      
      VolumeInfo worldVInfo = constructWorldVolume(_config);

      if ( _verbosityLevel > 0) {
          cout << __func__ << " worldVInfo.centerInParent : " <<  worldVInfo.centerInParent << endl;
          cout << __func__ << " worldVInfo.centerInWorld  : " <<  worldVInfo.centerInWorld  << endl;
      }

      VolumeInfo hallInfo  = constructHall(worldVInfo, _config);

      if ( _verbosityLevel > 0) {
          cout << __func__ << " hallInfo.centerInParent   : " <<  hallInfo.centerInParent << endl;
          cout << __func__ << " hallInfo.centerInWorld    : " <<  hallInfo.centerInWorld  << endl;
          cout << __func__ << " hallInfo.centerInMu2e()   : " <<  hallInfo.centerInMu2e() << endl;
      }
      
      // ************************  END: NO CHANGES NEEDED   ************************
      
      

      //CALLS constructExtMonFNAL, which has 4 cases of SetSensitiveDetector but no AddNew
      //SetSensitiveDet is called in loop, so need a vector of logical volumes
      constructProtonBeamDump(hallInfo, _config);//commented out SetSensitive, NEED TO FIX!!!!

      constructDS(hallInfo, _config);//NO CHANGES NEEDED
      
      
      //Here's my test case.  instead of being a void function, it returns a G4LogicalVolume*,
      //which I made a member function of this class.  Then, I access this member function
      //in constructSDandField() call to instantiateSensitiveDetectors in order to set the SensitiveDetectors
      psVacuumLogical_ = constructPS(hallInfo, _config);//CALLS SetSensitiveDetector 1x
    
  
      
      
        // ************************  BEGIN: NO CHANGES NEEDED   ************************
      
      constructPSEnclosure(hallInfo, _config);
      constructTS(hallInfo, _config);
      
        // ************************  BEGIN: NO CHANGES NEEDED   END? ************************

      
      //calls ConstructTTrackerv3 which has 4 cases of SetSensitiveDet, no AddNew
      //commented out SetSensitiveDet in ConstructTTrackerv3.cc, NEED TO FIX!!!!!
      //commented out SetSensitiveDet in ConstructTTrackerTDR.cc, NEED TO FIX!!!!!
    VolumeInfo trackerInfo = constructTracker();
      

      
      //calls constructStoppingTarget, which has 2 cases of SetSensitiveDet, no AddNew
      //I implemented this one below in instantiateSensitiveDetectors()
    VolumeInfo targetInfo  = constructTarget();
    
      
      //No AddNewDet, but
      //6 cases of SetSensitiveDetector:
      //4 in these types of structures: if (pabs->isAvailable(0))... inside //**** Modifiable MECO-Style proton absorber ****//
      //2 in these tyoes of structuresL if(paSD) protonabs1Info.logical->SetSensitiveDetector (paSD);, inside //**** MECO-Style proton absorber ****//
      //also, this file calls the file HelicalProtonAbsorber, which has 1 case of SetSensitiveDetector
    constructProtonAbsorber(_config);
      //commented out SetSensitiveDet in constructProtonAbsorber.cc & HelicalProtonAbsorber.cc, NEED TO FIX BOTH!!!!!

      
      //calls constructDiskCalorimater, which has 4 cases of SetSenstiveDet, no AddNew
      //commented out SetSensitiveDet in constructDiskCalorimeter.cc, NEED TO FIX!!!!!
    VolumeInfo calorimeterInfo = constructCal();

      
      
      
      // ************************  BEGIN: NO CHANGES NEEDED   ************************
      
    // This is just placeholder for now - and might be misnamed.
    constructMagnetYoke();

    // Check for stale names

    if ( _config.getBool("hasExternalShielding",false) ) {
      constructExternalShielding(hallInfo, _config);
    }

    // This is for pipes, cable runs, Electronics racks, etc.
    constructServicesGeom(hallInfo, _config);

    if ( _config.getBool("hasTSdA",false) ) {
      constructTSdA(_config);
    }

    if ( _config.getBool("hasMBS",false) ) {
      constructMBS(_config);
    }

      // ************************  END: NO CHANGES NEEDED   ************************
      
      
      
      
    if ( _config.getBool("mstm.build", false) ) {
      constructMSTM(hallInfo, _config);//1 case of SetSensitiveDet for SensitiveDetectorName::STMDet, no AddNew
        //commented out SetSensitiveDet in constructMSTM.cc, NEED TO FIX!!!!!
    }

    if ( _config.getBool("hasSTM",false) ) {
      constructSTM(_config);//2 cases of SetSensitiveDet like: if(sd1) stmDet1.logical->SetSensitiveDetector(sd1); no AddNew
        //commented out SetSensitiveDet in constructSTM.cc, NEED TO FIX!!!!!
    }
      
      
      // _geom is member data of Mu2eG4Universe, from which this inherits
      // it is a ref to a const GeometryService object
    if (  const_cast<GeometryService&>(_geom).hasElement<CosmicRayShield>() ) {

      //GeomHandle<CosmicRayShield> CosmicRayShieldGeomHandle;
      constructCRV(hallInfo,_config);//1 case of SetSensDet for SensitiveDetectorName::CRSScintillatorBar(), no AddNew
        //commented out SetSensitiveDet, NEED TO FIX!!!!!
    }

      
    /*

    ********  //44 INSTANCES!  TAKE THIS OUT FOR NOW.  IT RUNS WITHOUT THIS in place   ********
    constructVirtualDetectors(_config); // beware of the placement order of this function
     
     */

    constructVisualizationRegions(worldVInfo, _config);//STAYS, NO CHANGES

      //********  TAKE THIS OUT FOR NOW.  IT RUNS WITHOUT THIS in place   ********
      //this goes through a whole list of SDs and Adds them and Sets them - it will be an issue
      //it gets the list from the "sensitiveVolumes" list in the fcl file
      //for my configuration, this is empty
      //I will comment it out for now
      //maybe move to constructSDandField()?
    //sdHelper_->instantiateLVSDs(_config);
      

    if ( _verbosityLevel > 0) {
      mf::LogInfo log("GEOM");
      log << "Mu2e Origin:          " << worldGeom->mu2eOriginInWorld() << "\n";
    }

    //for MT, moved to ConstructSDAndField()
    // Create magnetic fields and managers only after all volumes have been defined.
    //constructBFieldAndManagers();
      
    constructStepLimiters();

    // Write out mu2e geometry into a gdml file.
    if (writeGDML_) {
      G4GDMLParser parser;
      parser.Write(gdmlFileName_, worldVInfo.logical);
    }

      //return worldVInfo.physical;
      worldVolume_ = worldVInfo.physical;
      return worldVolume_;
    
  }

  // Choose the selected tracker and build it.
  VolumeInfo Mu2eWorld::constructTracker(){

    // The tracker is built inside this volume.
    VolumeInfo const & detSolDownstreamVacInfo = _helper->locateVolInfo("DS3Vacuum");

    // z Position of the center of the DS solenoid parts, given in the Mu2e coordinate system.
    double z0DSdown = detSolDownstreamVacInfo.centerInMu2e().z();

    // Construct one of the trackers.
    VolumeInfo trackerInfo;

    if ( _config.getBool("hasTTracker",false) ) {
      int ver = _config.getInt("TTrackerVersion",3);
      if ( ver == 3 ){
        trackerInfo = constructTTrackerv3( detSolDownstreamVacInfo, _config );
      }
    } else {
      trackerInfo = constructDummyTracker( detSolDownstreamVacInfo.logical, z0DSdown, _config );
    }

    if ( _verbosityLevel > 0) {
      cout << "detSolDownstreamVacInfo.centerInMu2e().x() =" << detSolDownstreamVacInfo.centerInMu2e().x() << endl;
      cout << "detSolDownstreamVacInfo.centerInMu2e().y() =" << detSolDownstreamVacInfo.centerInMu2e().y() << endl;
      cout << "detSolDownstreamVacInfo.centerInMu2e().z() =" << detSolDownstreamVacInfo.centerInMu2e().z() << endl;
    }

    return trackerInfo;

  } // end Mu2eWorld::constructTracker

    
  // Either build the stopping target or a placeholder.
  // This should be call constrcutStoppingTarget but the name is already taken.
  VolumeInfo Mu2eWorld::constructTarget(){

    // The target is built inside this volume.
    VolumeInfo const & detSolUpstreamVacInfo   = ( _config.getBool("isDumbbell",false) ) ? _helper->locateVolInfo("DS3Vacuum") : _helper->locateVolInfo("DS2Vacuum");//DS3Vacuum to move the targets

    if ( _verbosityLevel > 0) {
      cout << "detSolUpstreamVacInfo.centerInWorld.z()=" << detSolUpstreamVacInfo.centerInWorld.z() << endl;
      cout << "detSolUpstreamVacInfo.centerInMu2e().z() =" << detSolUpstreamVacInfo.centerInMu2e().z() << endl;
    }

    // Buid the stopping target
    VolumeInfo targetInfo = ( _config.getBool("hasTarget",false) ) ?

      constructStoppingTarget( detSolUpstreamVacInfo,
                               _config )
      :

      constructDummyStoppingTarget( detSolUpstreamVacInfo,
                                    _config );
    return targetInfo;

  } // end Mu2eWorld::constructTarget

    
  // Construct the magnetic field managers and attach them to
  // the relevant volumes or to the world
  void Mu2eWorld::constructBFieldAndManagers(){
      
      if (G4Threading::G4GetThreadId() <= 0) {
          cout << "We are in Mu2eWorld::constructBFieldAndManagers() in Thread# "
               << G4Threading::G4GetThreadId() << endl;
      }

    GeomHandle<WorldG4> worldGeom;

    // Get some information needed further

    VolumeInfo const & ds2VacuumVacInfo = _helper->locateVolInfo("DS2Vacuum");
    VolumeInfo const & ds3VacuumVacInfo = _helper->locateVolInfo("DS3Vacuum");
    G4LogicalVolume* ds2Vacuum = ds2VacuumVacInfo.logical;
    G4LogicalVolume* ds3Vacuum = ds3VacuumVacInfo.logical;

    double ds2HL = static_cast<G4Tubs*>(ds2VacuumVacInfo.solid)->GetZHalfLength();
    G4ThreeVector ds2Z0 = ds2VacuumVacInfo.centerInWorld;
    G4ThreeVector beamZ0( ds2Z0.x(), ds2Z0.y(), ds2Z0.z()+ds2HL );

    // Decide on the G4 Stepper
    
      
      GeomHandle<BFieldConfig> bfConfig;
      
    bool needDSUniform = (bfConfig->dsFieldForm() == BFieldConfig::dsModelSplit || bfConfig->dsFieldForm() == BFieldConfig::dsModelUniform );
    bool needDSGradient = false;

    // Create field manager for the uniform DS field.
    if (needDSUniform) {

      _dsUniform = FieldMgr::forUniformField( bfConfig->getDSUniformValue()*CLHEP::tesla, worldGeom->mu2eOriginInWorld() );

      // Create field manager for the gradient field in DS3
      if(bfConfig->dsFieldForm() == BFieldConfig::dsModelSplit) {
        needDSGradient = true;
        _dsGradient = FieldMgr::forGradientField( bfConfig->getDSUniformValue().z()*CLHEP::tesla,
                                                  bfConfig->getDSGradientValue().z()*CLHEP::tesla/CLHEP::m,
                                                  beamZ0 );
      }
    }

    // Create global field managers; don't use FieldMgr here to avoid problem with ownership

    G4MagneticField * _field = new Mu2eGlobalField(worldGeom->mu2eOriginInWorld());
    G4Mag_EqRhs * _rhs  = new G4Mag_UsualEqRhs(_field);
    G4MagIntegratorStepper * _stepper;
    if ( g4stepperName_  == "G4ClassicalRK4" ) {
      _stepper = new G4ClassicalRK4(_rhs);
    } else if ( g4stepperName_  == "G4ClassicalRK4WSpin" ) {
      delete _rhs; // FIXME: avoid the delete
      _rhs  = new G4Mag_SpinEqRhs(_field);
      _stepper = new G4ClassicalRK4(_rhs, 12);
      if ( _verbosityLevel > 0) {
        cout << __func__ << " Replaced G4Mag_UsualEqRhs with G4ClassicalRK4WSpin "
             << "and used G4ClassicalRK4 with Spin" << endl;
      }
    } else if ( g4stepperName_  == "G4ImplicitEuler" ) {
      _stepper = new G4ImplicitEuler(_rhs);
    } else if ( g4stepperName_  == "G4ExplicitEuler" ) {
      _stepper = new G4ExplicitEuler(_rhs);
    } else if ( g4stepperName_  == "G4SimpleHeum" ) {
      _stepper = new G4SimpleHeum(_rhs);
    } else if ( g4stepperName_  == "G4HelixImplicitEuler" ) {
      _stepper = new G4HelixImplicitEuler(_rhs);
    } else if ( g4stepperName_  == "G4HelixSimpleRunge" ) {
      _stepper = new G4HelixSimpleRunge(_rhs);
    } else {
      _stepper = new G4SimpleRunge(_rhs);
    }
    G4ChordFinder * _chordFinder = new G4ChordFinder(_field,1.0e-2*CLHEP::mm,_stepper);
    G4FieldManager * _manager = new G4FieldManager(_field,_chordFinder,true);

    // G4TransportationManager takes ownership of _manager
    G4TransportationManager::GetTransportationManager()->SetFieldManager(_manager);

    // Define uniform field region in the detector solenoid, if neccessary
    if (bfConfig->dsFieldForm() == BFieldConfig::dsModelUniform  ){
      ds2Vacuum->SetFieldManager( _dsUniform->manager(), true);
      if ( _verbosityLevel > 0 ) cout << "Use uniform field in DS2" << endl;
    }
    if (bfConfig->dsFieldForm() == BFieldConfig::dsModelUniform || bfConfig->dsFieldForm() == BFieldConfig::dsModelSplit){
      if( needDSGradient ) {
        ds3Vacuum->SetFieldManager( _dsGradient->manager(), true);
      if ( _verbosityLevel > 0 ) cout << "Use gradient field in DS3" << endl;
      } else {
        ds3Vacuum->SetFieldManager( _dsUniform->manager(), true);
        if ( _verbosityLevel > 0 ) cout << "Use uniform field in DS3" << endl;
      }

      // if( _config.getBool("hasMBS",false) ) {
      //   VolumeInfo const & MBSMotherInfo = _helper->locateVolInfo("MBSMother");
      //   G4LogicalVolume* mbsMother = MBSMotherInfo.logical;
      //   mbsMother->SetFieldManager(0,true);
      // }

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


  namespace {

    // A helper function for Mu2eWorld::constructStepLimiters().
    // Find all logical volumes matching a wildcarded name and add steplimiters to them.
    void stepLimiterHelper ( std::string const& regexp, G4UserLimits* stepLimit ){
      boost::regex expression(regexp.c_str());
      std::vector<mu2e::VolumeInfo const*> vols =
        art::ServiceHandle<mu2e::G4Helper>()->locateVolInfo(expression);
      for ( auto v : vols ){
        v->logical->SetUserLimits( stepLimit );
      }
    }
  }

  // Adding a step limiter is a two step process.
  // 1) In the physics list constructor add a G4StepLimiter to the list of discrete
  //    physics processes attached to each particle species of interest.
  //
  // 2) In this code, create a G4UserLimits object and attach it to the logical
  //    volumes of interest.
  // The net result is specifying a step limiter for pairs of (logical volume, particle species).
  //
  void Mu2eWorld::constructStepLimiters(){

    G4LogicalVolume* ds2Vacuum      = _helper->locateVolInfo("DS2Vacuum").logical;
    G4LogicalVolume* ds3Vacuum      = _helper->locateVolInfo("DS3Vacuum").logical;
    G4LogicalVolume* tracker        = _helper->locateVolInfo("TrackerMother").logical;
    G4LogicalVolume* caloMother     = _helper->locateVolInfo("CalorimeterMother").logical;
    G4LogicalVolume* stoppingtarget = _helper->locateVolInfo("StoppingTargetMother").logical;

    vector<G4LogicalVolume*> psVacua;
    psVacua.push_back( _helper->locateVolInfo("PSVacuum").logical );

    vector<G4LogicalVolume*> tsVacua;
    tsVacua.push_back( _helper->locateVolInfo("TS1Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS2Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS3Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS4Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS5Vacuum").logical );

    vector<G4LogicalVolume*> mbsLVS;
    mbsLVS.push_back( _helper->locateVolInfo("MBSMother").logical );

    // We may make separate G4UserLimits objects per logical volume but we choose not to.
    // At some it might be interesting to make several step limiters, each with different
    // limits.  For now that is not necessary.
    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
    G4UserLimits* stepLimit = reg.add( G4UserLimits(bfieldMaxStep_) );

    // Add the step limiters to the interesting volumes.
    // Keep them separated so that we can add different step limits should we decide to.
    ds2Vacuum->SetUserLimits( stepLimit );
    ds3Vacuum->SetUserLimits( stepLimit );

    tracker->SetUserLimits( stepLimit );
    caloMother->SetUserLimits( stepLimit );
    stoppingtarget->SetUserLimits( stepLimit );

    for ( auto lv : psVacua ){
      lv->SetUserLimits( stepLimit);
    }

    for ( auto lv : tsVacua ){
      lv->SetUserLimits( stepLimit);
    }

    for ( auto lv : mbsLVS ){
      lv->SetUserLimits( stepLimit );
    }

    // Now do all of the tracker related envelope volumes, using regex's with wildcards.
    stepLimiterHelper("^TTrackerPlaneEnvelope_.*$",                 stepLimit );
    stepLimiterHelper("^TTrackerSupportServiceEnvelope_.*$",        stepLimit );
    stepLimiterHelper("^TTrackerSupportServiceSectionEnvelope_.*$", stepLimit );

    // An option to limit the step size in these non-vaccum volumes to
    // visually validate geometry of the filter channel
    //
    // ==== BEGIN COMMENT-OUT: to allow construction of new building dirt volumes w/o overlaps (knoepfel)
    //
    //     G4LogicalVolume* emfcMagnetAperture = _helper->locateVolInfo("ExtMonFNALfilterMagnetAperture").logical;
    //     if(emfcMagnetAperture) {
    //       const double maxStepLength = _config.getDouble("extMonFNAL.maxG4StepLength", 0)*CLHEP::millimeter;
    //       if(maxStepLength > 0) {
    //         std::cout<<"Adding step limiter for ExtMonFNALFilterMagnet: maxStepLength = "<<maxStepLength<<std::endl;
    //         G4UserLimits* emfcStepLimit = reg.add(G4UserLimits(maxStepLength));
    //         emfcMagnetAperture->SetUserLimits(emfcStepLimit);
    //         _helper->locateVolInfo("ExtMonFNALfilterMagnetIron").logical->SetUserLimits(emfcStepLimit);
    //       }
    //     }
    //
    // ==== END COMMENT-OUT

  } // end Mu2eWorld::constructStepLimiters(){


  // Construct calorimeter if needed.
  VolumeInfo Mu2eWorld::constructCal(){

        // The calorimeter is built inside this volume.
        VolumeInfo const & detSolDownstreamVacInfo = _helper->locateVolInfo("DS3Vacuum");

        // Construct one of the calorimeters.
        VolumeInfo calorimeterInfo;
        if ( _config.getBool("hasDiskCalorimeter",false) ) {
           calorimeterInfo = constructDiskCalorimeter( detSolDownstreamVacInfo,_config );
        }

    return calorimeterInfo;

  } // end Mu2eWorld::constructCal


  // A place holder for now.
  void Mu2eWorld::constructMagnetYoke(){
  } // end Mu2eWorld::constructMagnetYoke

    
  // instantiateSensitiveDetectors
  void Mu2eWorld::instantiateSensitiveDetectors(){

      //art::ServiceHandle<GeometryService> geom;//I don't see this being used!
      G4SDManager* SDman = G4SDManager::GetSDMpointer();
      
      //get the LVStore singleton
      G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();


      if (G4Threading::G4GetThreadId() <= 0){
          
          std::cout << "in Thread" << G4Threading::G4GetThreadId() << ", Mu2eWorld::instantiateSDs() data products to be produced: " << std::endl;
          vector<string> const& instanceNamesTEST = sdHelper_->stepInstanceNamesToBeProduced();
     
          for (unsigned i = 0 ; i<instanceNamesTEST.size(); i++) {
              cout << instanceNamesTEST.at(i) << endl;
          }
          
      }

      
      // G4 takes ownership and will delete the detectors at the job end
    
      //done
    if(sdHelper_->enabled(StepInstanceName::tracker)) {
        
        //cout << "in Mu2eWorld::instantiateSDs(), tracker is enabled. setting StrawSDs" << std::endl;
        
        StrawSD* strawSD =
        new StrawSD( SensitiveDetectorName::TrackerGas(), _config );
        SDman->AddNewDetector(strawSD);
        
        //std::cout << "for SD: " << SensitiveDetectorName::TrackerGas() << std::endl;

        //loop over all of the LV names and find the ones we need
        //set the SensitiveDetectors for these
        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
            
            G4String LVname = (*pos)->GetName();
            
            if (LVname.find("TTrackerStrawGas_") != std::string::npos  ) {
                //if ((*pos)->GetName().find(basename) != std::string::npos  ) {
                //std::cout << "Setting SD for LV: " << (*pos)->GetName() << std::endl;
                store->GetVolume(LVname)->SetSensitiveDetector(strawSD);
            }
        }

     }//if tracker

     //done
    if(sdHelper_->enabled(StepInstanceName::ttrackerDS)) {
        
        TTrackerPlaneSupportSD* ttdsSD =
        new TTrackerPlaneSupportSD( SensitiveDetectorName::TTrackerPlaneSupport(), _config );
        SDman->AddNewDetector(ttdsSD);
        
        //std::cout << "for SD: " << SensitiveDetectorName::TTrackerPlaneSupport() << std::endl;
        
        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
            
            G4String LVname = (*pos)->GetName();
            
            if (LVname.find("TTrackerPlaneSupport_") != std::string::npos  ) {
                //std::cout << "Setting SD for LV: " << (*pos)->GetName() << std::endl;
                store->GetVolume(LVname)->SetSensitiveDetector(ttdsSD);
            }
        }
        
    }//if ttrackerDS
      
      
      
    if(sdHelper_->enabled(StepInstanceName::virtualdetector)) {
        
        Mu2eSensitiveDetector* vdSD =
            new Mu2eSensitiveDetector( SensitiveDetectorName::VirtualDetector(), _config );
        SDman->AddNewDetector(vdSD);

      if (activeWr_Wl_SD_) {//see constructTTrackerv3.cc
          
          //done
        if (sdHelper_->enabled(StepInstanceName::trackerSWires)) {
            TrackerWireSD *ttwsSD =
                new TrackerWireSD(SensitiveDetectorName::TrackerSWires(), _config);
            SDman->AddNewDetector(ttwsSD);
            
            //trackerSWires need to be fixed if virtual detector is fixed
            //std::cout << "for SD: " << SensitiveDetectorName::TrackerSWires() << std::endl;

            for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
                
                G4String LVname = (*pos)->GetName();
                
                if (LVname.find("TTrackerStrawWire_") != std::string::npos  ) {
                    //std::cout << "Setting SD for LV: " << (*pos)->GetName() << std::endl;
                    store->GetVolume(LVname)->SetSensitiveDetector(ttwsSD);
                }
            }
        }
          
          //done
        if (sdHelper_->enabled(StepInstanceName::trackerWalls)) {//see constructTTrackerv3.cc
            TrackerWireSD *ttwlSD =
                new TrackerWireSD(SensitiveDetectorName::TrackerWalls(), _config);
            SDman->AddNewDetector(ttwlSD);
            //trackerWalls need to be fixed if virtual detector is fixed
            
            //std::cout << "for SD: " << SensitiveDetectorName::TrackerWalls() << std::endl;
            for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
                
                G4String LVname = (*pos)->GetName();
                
                if (LVname.find("TTrackerStrawWall_") != std::string::npos  ) {
                    //std::cout << "Setting SD for LV: " << (*pos)->GetName() << std::endl;
                    store->GetVolume(LVname)->SetSensitiveDetector(ttwlSD);
                }
            }
        }
          
      }//if activeWr_Wl_SD_

    }//if virtualdetector

      
      
      // _geom is member data of Mu2eG4Universe, from which this inherits (GeoService const&)
      //done
    if ( const_cast<GeometryService&>(_geom).hasElement<Calorimeter>() ) {
        
      if(sdHelper_->enabled(StepInstanceName::calorimeter)) {
          
          CaloCrystalSD* ccSD =
            new CaloCrystalSD( SensitiveDetectorName::CaloCrystal(), _config );
          SDman->AddNewDetector(ccSD);
          
          //std::cout << "for SD: " << SensitiveDetectorName::CaloCrystal() << std::endl;
          
          for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
              
              G4String LVname = (*pos)->GetName();
              
              if (LVname.find("CrystalLog") != std::string::npos  ) {
                  //std::cout << "Setting SD for LV: " << (*pos)->GetName() << std::endl;
                  store->GetVolume(LVname)->SetSensitiveDetector(ccSD);
              }
          }
      }

      if(sdHelper_->enabled(StepInstanceName::calorimeterRO)) {
          CaloReadoutSD* crSD =
            new CaloReadoutSD( SensitiveDetectorName::CaloReadout(), _config );
          SDman->AddNewDetector(crSD);
        
          //std::cout << "for SD: " << SensitiveDetectorName::CaloReadout() << std::endl;
          
          for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
              
              G4String LVname = (*pos)->GetName();
              
              if (LVname.find("CrystalROLog") != std::string::npos  ) {
                  //std::cout << "Setting SD for LV: " << (*pos)->GetName() << std::endl;
                  store->GetVolume(LVname)->SetSensitiveDetector(crSD);
              }
          }
      }

      if(sdHelper_->enabled(StepInstanceName::calorimeterROCard)) {
          CaloReadoutCardSD* crCardSD =
            new CaloReadoutCardSD( SensitiveDetectorName::CaloReadoutCard(), _config );
          SDman->AddNewDetector(crCardSD);
          
          //std::cout << "for SD: " << SensitiveDetectorName::CaloReadoutCard() << std::endl;
          
          for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
              
              G4String LVname = (*pos)->GetName();
              
              if (LVname.find("ElectronicsROLog") != std::string::npos  ) {
                  //std::cout << "Setting SD for LV: " << (*pos)->GetName() << std::endl;
                  store->GetVolume(LVname)->SetSensitiveDetector(crCardSD);
              }
          }
      }
      
      if(sdHelper_->enabled(StepInstanceName::calorimeterCrate)) {
          CaloCrateSD* cCrateSD =
            new CaloCrateSD( SensitiveDetectorName::CaloCrate(), _config );
          SDman->AddNewDetector(cCrateSD);
          
          
          //std::cout << "for SD: " << SensitiveDetectorName::CaloCrate() << std::endl;
          
          for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
              
              G4String LVname = (*pos)->GetName();
              
              if (LVname.find("activeStripBoardLog") != std::string::npos  ) {
                  //std::cout << "Setting SD for LV: " << (*pos)->GetName() << std::endl;
                  store->GetVolume(LVname)->SetSensitiveDetector(cCrateSD);
              }
          }
      }

    }//if calorimeter

    if(true) { // this SD does not derive from Mu2eSensitiveDetector as it does not produce StepPointMCCollection
        GeomHandle<mu2e::ExtMonFNAL::ExtMon> extmon;
        SDman->AddNewDetector(new ExtMonFNALPixelSD(_config, *extmon));
    }

      //done
    if(sdHelper_->enabled(StepInstanceName::stoppingtarget)) {
        
        Mu2eSensitiveDetector* stSD =
            new Mu2eSensitiveDetector( SensitiveDetectorName::StoppingTarget(), _config );
        SDman->AddNewDetector(stSD);
        
        //loop over all of the LV names to find ones we need
        //set the SensitiveDetectors for these
        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
            
            G4String LVname = (*pos)->GetName();
           
            if (LVname.find("Foil_") != std::string::npos  ) {
                //std::cout << "Setting SD for LV: " << (*pos)->GetName() << std::endl;
                store->GetVolume(LVname)->SetSensitiveDetector(stSD);
            }
            
            if (LVname.find("FoilSupportStructure_") != std::string::npos  ) {
                //std::cout << "Setting SD for LV: " << (*pos)->GetName() << std::endl;
                store->GetVolume(LVname)->SetSensitiveDetector(stSD);
            }
        }
        
    }//if stoppingtarget

      
      
    if(sdHelper_->enabled(StepInstanceName::CRV)) {
        
      Mu2eSensitiveDetector* sbSD =
        new Mu2eSensitiveDetector( SensitiveDetectorName::CRSScintillatorBar(), _config);
      SDman->AddNewDetector(sbSD);
    }

    if(sdHelper_->enabled(StepInstanceName::protonabsorber)) {
      Mu2eSensitiveDetector* paSD =
        new Mu2eSensitiveDetector(    SensitiveDetectorName::ProtonAbsorber(),  _config);
      SDman->AddNewDetector(paSD);
    }

      //done
    if(sdHelper_->enabled(StepInstanceName::PSVacuum)) {
        
        Mu2eSensitiveDetector* psVacuumSD =
            new Mu2eSensitiveDetector( SensitiveDetectorName::PSVacuum(), _config );
      
        SDman->AddNewDetector(psVacuumSD);
        
        bool psVacuumSensitive = _config.getBool("PS.Vacuum.Sensitive", false);
        
        if(psVacuumSensitive) {
            psVacuumLogical_->SetSensitiveDetector(psVacuumSD);
        }
    }
      

    if(sdHelper_->enabled(StepInstanceName::STMDet)) {
      Mu2eSensitiveDetector* STMDetSD =
        new Mu2eSensitiveDetector(    SensitiveDetectorName::STMDet(),  _config);
      SDman->AddNewDetector(STMDetSD);
    }

    if(sdHelper_->enabled(StepInstanceName::panelEBKey)) {
      Mu2eSensitiveDetector* EBKeySD =
        new Mu2eSensitiveDetector(    SensitiveDetectorName::panelEBKey(),  _config);
      SDman->AddNewDetector(EBKeySD);
    }

  } // instantiateSensitiveDetectors

} // end namespace mu2e
