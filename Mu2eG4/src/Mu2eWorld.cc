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
#include "Mu2eG4/inc/constructSaddles.hh"
#include "Mu2eG4/inc/constructServicesGeom.hh"
#include "Mu2eG4/inc/constructTSdA.hh"
#include "Mu2eG4/inc/constructMBS.hh"
#include "Mu2eG4/inc/constructVirtualDetectors.hh"
#include "Mu2eG4/inc/constructVirtualDetectorSDs.hh"
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
#include "Mu2eG4/inc/TrackerPlaneSupportSD.hh"
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
#include "Mu2eG4/inc/constructTracker.hh"
#include "Mu2eG4/inc/constructStoppingTarget.hh"
#include "Mu2eG4/inc/constructDummyStoppingTarget.hh"
#include "Mu2eG4/inc/constructDiskCalorimeter.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "TrackerGeom/inc/Tracker.hh"
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
#include "G4PropagatorInField.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4UserLimits.hh"
#include "G4ClassicalRK4.hh"
#include "G4ImplicitEuler.hh"
#include "G4ExplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#if G4VERSION>4103
#include "G4DormandPrince745.hh"
#include "G4BogackiShampine23.hh"
#endif
#include "G4GDMLParser.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"

#include "Mu2eG4/inc/Mu2eGlobalField.hh"

#include "boost/regex.hpp"

using namespace std;

namespace mu2e {

  Mu2eWorld::Mu2eWorld(SensitiveDetectorHelper *sdHelper/*no ownership passing*/)
    : sdHelper_(sdHelper)
    , pset_(fhicl::ParameterSet())
    , activeWr_Wl_SD_(_config.getBool("tracker.ActiveWr_Wl_SD",false))
    , writeGDML_(_config.getBool("writeGDML",false))
    , gdmlFileName_(_config.getString("GDMLFileName","mu2e.gdml"))
    , g4stepperName_(_config.getString("g4.stepper","G4SimpleRunge"))
    , g4epsilonMin_(_config.getDouble("g4.epsilonMin"))
    , g4epsilonMax_(_config.getDouble("g4.epsilonMax"))
    , g4DeltaOneStep_(_config.getDouble("g4.deltaOneStep")*CLHEP::mm)
    , g4DeltaIntersection_(_config.getDouble("g4.deltaIntersection")*CLHEP::mm)
    , g4DeltaChord_(_config.getDouble("g4.deltaChord")*CLHEP::mm)
    , bfieldMaxStep_(_config.getDouble("bfield.maxStep", 20.)*CLHEP::mm)
    , strawGasMaxStep_(_config.getDouble("strawGas.maxStep", 20.)*CLHEP::mm)
    , useEmOption4InTracker_(_config.getBool("g4.useEmOption4InTracker",false))
  {
    _verbosityLevel = _config.getInt("world.verbosityLevel", 0);
    _g4VerbosityLevel = _config.getInt("g4.diagLevel", 0);
  }

  Mu2eWorld::Mu2eWorld(const fhicl::ParameterSet& pset,
                       SensitiveDetectorHelper *sdHelper/*no ownership passing*/)
    : sdHelper_(sdHelper)
    , pset_(pset)
    , activeWr_Wl_SD_(true)
    , writeGDML_(pset.get<bool>("debug.writeGDML"))
    , gdmlFileName_(pset.get<std::string>("debug.GDMLFileName"))
    , g4stepperName_(pset.get<std::string>("physics.stepper"))
    , g4epsilonMin_(pset.get<double>("physics.epsilonMin"))
    , g4epsilonMax_(pset.get<double>("physics.epsilonMax"))
    , g4DeltaOneStep_(pset.get<double>("physics.deltaOneStep")*CLHEP::mm)
    , g4DeltaIntersection_(pset.get<double>("physics.deltaIntersection")*CLHEP::mm)
    , g4DeltaChord_(pset.get<double>("physics.deltaChord")*CLHEP::mm)
    , g4StepMinimum_(pset.get<double>("physics.stepMinimum")*CLHEP::mm)
    , g4MaxIntSteps_(pset.get<int>("physics.maxIntSteps"))
    , bfieldMaxStep_(pset.get<double>("physics.bfieldMaxStep")*CLHEP::mm)
    , strawGasMaxStep_(pset.get<double>("physics.strawGasMaxStep")*CLHEP::mm)
    , limitStepInAllVolumes_(pset.get<bool>("physics.limitStepInAllVolumes"))
    , useEmOption4InTracker_(pset.get<bool>("physics.useEmOption4InTracker",false))
  {
    _verbosityLevel = pset.get<int>("debug.worldVerbosityLevel");
    _g4VerbosityLevel = pset.get<int>("debug.diagLevel");
  }


    // This is the callback called by G4 via G4VPhysicalVolume* WorldMaker::Construct()
  G4VPhysicalVolume * Mu2eWorld::construct(){

        // Construct all of the Mu2e world, hall, detectors, beamline ...
        return constructWorld();
  }

    // This is the callback called by G4 via void WorldMaker::ConstructSDandField()
  void Mu2eWorld::constructSDandField(){

      sdHelper_->instantiateLVSDs(_config);
      instantiateSensitiveDetectors();
      constructBFieldAndManagers();

  }


  // Construct all of the Mu2e world, hall, detectors, beamline ...
  G4VPhysicalVolume * Mu2eWorld::constructWorld(){

      // If you play with the order of these calls, you may break things.
      GeomHandle<WorldG4> worldGeom;
      G4ThreeVector tmpTrackercenter = GeomHandle<DetectorSystem>()->getOrigin();

      if (activeWr_Wl_SD_) {
          TrackerWireSD::setMu2eDetCenterInWorld( tmpTrackercenter );
      }

      VolumeInfo worldVInfo = constructWorldVolume(_config);

      if ( _verbosityLevel > 0) {
          G4cout << __func__ << " worldVInfo.centerInParent : " <<  worldVInfo.centerInParent << G4endl;
          G4cout << __func__ << " worldVInfo.centerInWorld  : " <<  worldVInfo.centerInWorld  << G4endl;
      }

      VolumeInfo hallInfo  = constructHall(worldVInfo, _config);

      if ( _verbosityLevel > 0) {
          G4cout << __func__ << " hallInfo.centerInParent   : " <<  hallInfo.centerInParent << G4endl;
          G4cout << __func__ << " hallInfo.centerInWorld    : " <<  hallInfo.centerInWorld  << G4endl;
          G4cout << __func__ << " hallInfo.centerInMu2e()   : " <<  hallInfo.centerInMu2e() << G4endl;
      }

      constructProtonBeamDump(hallInfo, _config);
      constructDS(hallInfo, _config);

      //Here's a test case for a another way of making SDs active.
      //Instead of being a void function, it returns a G4LogicalVolume*,
      //which I made a member function of this class.  Then, I access this member function
      //in constructSDandField() call to instantiateSensitiveDetectors in order to set the SensitiveDetectors
      psVacuumLogical_ = constructPS(hallInfo, _config);

      constructPSEnclosure(hallInfo, _config);
      constructTS(hallInfo, _config);
      VolumeInfo trackerInfo = constructTracker();
      VolumeInfo targetInfo  = constructTarget();
      constructProtonAbsorber(_config);
      VolumeInfo calorimeterInfo = constructCal();

      // This is just placeholder for now - and might be misnamed.
      constructMagnetYoke();

      if ( _config.getBool("hasExternalShielding",false) ) {
          constructExternalShielding(hallInfo, _config);
      }

      // This is for saddles holding up cryostats
      constructSaddles(hallInfo, _config);

      // This is for pipes, cable runs, Electronics racks, etc.
      constructServicesGeom(hallInfo, _config);

      if ( _config.getBool("hasTSdA",false) ) {
          constructTSdA(_config);
      }

      if ( _config.getBool("hasMBS",false) ) {
          constructMBS(_config);
      }

      if ( _config.getBool("mstm.build", false) ) {
          constructMSTM(hallInfo, _config);
      }

      if ( _config.getBool("hasSTM",false) ) {
          constructSTM(_config);
      }

      // _geom is member data of Mu2eG4Universe, from which this inherits
      // it is a ref to a const GeometryService object
      if (  const_cast<GeometryService&>(_geom).hasElement<CosmicRayShield>() ) {
          GeomHandle<CosmicRayShield> CosmicRayShieldGeomHandle;
          constructCRV(hallInfo,_config);
      }

      constructVirtualDetectors(_config); // beware of the placement order of this function

      constructVisualizationRegions(worldVInfo, _config);

      if ( _verbosityLevel > 0) {
          mf::LogInfo log("GEOM");
          log << "Mu2e Origin:          " << worldGeom->mu2eOriginInWorld() << "\n";
      }

      // creating regions to be able to asign special cut and EM options

      const fhicl::ParameterSet& minRangeRegionCutsPSet{
        pset_.get<fhicl::ParameterSet>("physics.minRangeRegionCuts",fhicl::ParameterSet())};

      if (!minRangeRegionCutsPSet.is_empty()) {
        const std::vector<std::string> regionNames{minRangeRegionCutsPSet.get_names()};
        for(const auto& regionName : regionNames) {
          G4Region* region = new G4Region(regionName); // G4RegionStore takes ownership
          VolumeInfo const & volInfo = _helper->locateVolInfo(regionName);
          volInfo.logical->SetRegion(region);
          region->AddRootLogicalVolume(volInfo.logical);

          G4ProductionCuts* regionProductionCuts = new G4ProductionCuts();
          G4double productionCut = minRangeRegionCutsPSet.get<double>(regionName);
          regionProductionCuts->SetProductionCut(productionCut);
          // the above sets the same cut for gamma, e- and e+, proton/ions
          G4double protonProductionCut = pset_.get<double>("physics.protonProductionCut");
          regionProductionCuts->SetProductionCut(protonProductionCut,"proton");
          region->SetProductionCuts(regionProductionCuts);

          if ( _verbosityLevel > 0 ) {
            G4cout << __func__ << " Setting gamma, e- and e+ production cut for "
                      << regionName << " to " << productionCut << " mm and for proton to "
                      << protonProductionCut << " mm" << G4endl;
            G4cout << __func__ << " Resulting cuts for gamma, e-, e+, proton: ";
            for (auto const& rcut : regionProductionCuts->GetProductionCuts() ) {
              G4cout << " " << rcut;
            }
            G4cout << G4endl;
          }

        }
      }

      // special case for the tracker when we need a region to set a
      // different EM option even when no production cuts are set explicitly

      if ( useEmOption4InTracker_
           && !pset_.has_key("physics.minRangeRegionCuts.TrackerMother")) {
        G4Region* region = new G4Region("TrackerMother");
        trackerInfo.logical->SetRegion(region);
        region->AddRootLogicalVolume(trackerInfo.logical);
      }

      constructStepLimiters();

      // Write out mu2e geometry into a gdml file.
      if (writeGDML_) {
          G4GDMLParser parser;
          parser.Write(gdmlFileName_, worldVInfo.logical);
      }

      return worldVInfo.physical;

  }//Mu2eWorld::constructWorld()


  // Choose the selected tracker and build it.
  VolumeInfo Mu2eWorld::constructTracker(){

    // The tracker is built inside this volume.
    std::string theDS3("DS3Vacuum");
    if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";
    VolumeInfo const & detSolDownstreamVacInfo = _helper->locateVolInfo(theDS3);

    // Construct one of the trackers.
    VolumeInfo trackerInfo;

    if ( _config.getBool("hasTracker",false) ) {
      if ( _config.getInt("TrackerVersion",3)  != 5 ){
        throw cet::exception("GEOM")
          << "Current code only supports TrackerVersion = 5\n"
          << "Requested version: "
          << _config.getInt("TrackerVersion",3);
      }
      trackerInfo = constructTrackerv5( detSolDownstreamVacInfo, _config);
    }

    if ( _verbosityLevel > 0) {
      G4cout << __func__ << "detSolDownstreamVacInfo.centerInMu2e().x() =" << detSolDownstreamVacInfo.centerInMu2e().x() << G4endl;
      G4cout << __func__ << "detSolDownstreamVacInfo.centerInMu2e().y() =" << detSolDownstreamVacInfo.centerInMu2e().y() << G4endl;
      G4cout << __func__ << "detSolDownstreamVacInfo.centerInMu2e().z() =" << detSolDownstreamVacInfo.centerInMu2e().z() << G4endl;
    }

    return trackerInfo;

  }//Mu2eWorld::constructTracker


  // Either build the stopping target or a placeholder.
  // This should be call constrcutStoppingTarget but the name is already taken.
  VolumeInfo Mu2eWorld::constructTarget(){

    // The target is built inside this volume.
    VolumeInfo const & detSolUpstreamVacInfo   = ( _config.getBool("isDumbbell",false) ) ? _helper->locateVolInfo("DS3Vacuum") : _helper->locateVolInfo("DS2Vacuum");//DS3Vacuum to move the targets

    if ( _verbosityLevel > 0) {
      G4cout << __func__ << "detSolUpstreamVacInfo.centerInWorld.z()=" << detSolUpstreamVacInfo.centerInWorld.z() << G4endl;
      G4cout << __func__ << "detSolUpstreamVacInfo.centerInMu2e().z() =" << detSolUpstreamVacInfo.centerInMu2e().z() << G4endl;
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
    if ( _g4VerbosityLevel > 0 ) G4cout << __func__ << " Setting up " << g4stepperName_ << " stepper" << G4endl;
    if ( g4stepperName_  == "G4ClassicalRK4" ) {
      _stepper = new G4ClassicalRK4(_rhs);
    } else if ( g4stepperName_  == "G4ClassicalRK4WSpin" ) {
      delete _rhs; // FIXME: avoid the delete
      _rhs  = new G4Mag_SpinEqRhs(_field);
      _stepper = new G4ClassicalRK4(_rhs, 12);
      if ( _g4VerbosityLevel > 0) {
        G4cout << __func__ << " Replaced G4Mag_UsualEqRhs with G4ClassicalRK4WSpin "
             << "and used G4ClassicalRK4 with Spin" << G4endl;
      }
    } else if ( g4stepperName_  == "G4DormandPrince745WSpin" ) {
      delete _rhs; // FIXME: avoid the delete
      _rhs  = new G4Mag_SpinEqRhs(_field);
      _stepper = new G4DormandPrince745(_rhs, 12);
      if ( _g4VerbosityLevel > 0) {
        G4cout << __func__ << " Replaced G4Mag_UsualEqRhs with G4DormandPrince745WSpin "
             << "and used G4DormandPrince745 with Spin" << G4endl;
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
#if G4VERSION>4103
    } else if ( g4stepperName_  == "G4DormandPrince745" ) {
      _stepper = new G4DormandPrince745(_rhs);
    } else if ( g4stepperName_  == "G4BogackiShampine23" ) {
      _stepper = new G4BogackiShampine23(_rhs);
#endif
    } else {
      _stepper = new G4SimpleRunge(_rhs);
      if ( _g4VerbosityLevel > -1 ) G4cout << __func__ << " Using default G4SimpleRunge stepper" << G4endl;
    }
    G4ChordFinder * _chordFinder = new G4ChordFinder(_field,g4StepMinimum_,_stepper);
    G4FieldManager * _manager = new G4FieldManager(_field,_chordFinder,true);

    // G4TransportationManager takes ownership of _manager

    G4TransportationManager* transporationMgr =
      G4TransportationManager::GetTransportationManager();
    transporationMgr->SetFieldManager(_manager);

    // Define uniform field region in the detector solenoid, if neccessary
    if (bfConfig->dsFieldForm() == BFieldConfig::dsModelUniform  ){
      ds2Vacuum->SetFieldManager( _dsUniform->manager(), true);
      if ( _verbosityLevel > 0 ) G4cout << __func__ << " Use uniform field in DS2" << G4endl;
    }
    if (bfConfig->dsFieldForm() == BFieldConfig::dsModelUniform || bfConfig->dsFieldForm() == BFieldConfig::dsModelSplit){
      if( needDSGradient ) {
        ds3Vacuum->SetFieldManager( _dsGradient->manager(), true);
      if ( _verbosityLevel > 0 ) G4cout << __func__ << " Use gradient field in DS3" << G4endl;
      } else {
        ds3Vacuum->SetFieldManager( _dsUniform->manager(), true);
        if ( _verbosityLevel > 0 ) G4cout << __func__ << " Use uniform field in DS3" << G4endl;
      }

      // if( _config.getBool("hasMBS",false) ) {
      //   VolumeInfo const & MBSMotherInfo = _helper->locateVolInfo("MBSMother");
      //   G4LogicalVolume* mbsMother = MBSMotherInfo.logical;
      //   mbsMother->SetFieldManager(0,true);
      // }

    }

    // Adjust properties of the integrators to control accuracy vs time.

    if ( _dsUniform.get() != 0 ){
      G4double deltaIntersection = 0.00001*CLHEP::mm;
      _dsUniform->manager()->SetDeltaIntersection(deltaIntersection);
    }

    if ( _dsGradient.get() != 0 ){
      G4double deltaIntersection = 0.00001*CLHEP::mm;
      _dsGradient->manager()->SetDeltaIntersection(deltaIntersection);
    }

    _manager->SetMinimumEpsilonStep(g4epsilonMin_);
    _manager->SetMaximumEpsilonStep(g4epsilonMax_);
    _manager->SetDeltaOneStep(g4DeltaOneStep_);
    _manager->SetDeltaIntersection(g4DeltaIntersection_);
    _chordFinder->SetDeltaChord(g4DeltaChord_);

    G4PropagatorInField* _propInField = transporationMgr->GetPropagatorInField();
    _propInField->SetMaxLoopCount(g4MaxIntSteps_);


    if ( _g4VerbosityLevel > 0 ) {
      G4cout << __func__ << " Stepper precision parameters: " << G4endl;
      G4cout << __func__ << " g4epsilonMin        " << _manager->GetMinimumEpsilonStep() << G4endl;
      G4cout << __func__ << " g4epsilonMax        " << _manager->GetMaximumEpsilonStep() << G4endl;
      G4cout << __func__ << " g4DeltaOneStep      " << _manager->GetDeltaOneStep() << G4endl;
      G4cout << __func__ << " g4DeltaIntersection " << _manager->GetDeltaIntersection() << G4endl;
      G4cout << __func__ << " g4DeltaChord        " << _chordFinder->GetDeltaChord() << G4endl;
      G4cout << __func__ << " g4StepMinimum       "
      // << dynamic_cast<G4MagInt_Driver*>(_chordFinder->GetIntegrationDriver())->GetHmin() << G4endl;
      // the above assumes G4ChordFinder is instantiated in the way it is done above with the 3 parameters
      // does not work in 10.5+; fixme
	     << g4StepMinimum_ << G4endl;
      G4cout << __func__ << " g4MaxIntStep        " << _propInField->GetMaxLoopCount() << G4endl;
    }

  } // end Mu2eWorld::constructBFieldAndManagers


    // A helper function for Mu2eWorld::constructStepLimiters().
    // Find all logical volumes matching a wildcarded name and add steplimiters to them.
  void Mu2eWorld::stepLimiterHelper ( std::string const& regexp, G4UserLimits* stepLimit ) {
    boost::regex expression(regexp.c_str());
    std::vector<mu2e::VolumeInfo const*> vols =
      art::ServiceHandle<mu2e::G4Helper>()->locateVolInfo(expression);
    for ( auto v : vols ){
      v->logical->SetUserLimits( stepLimit );
      if(_g4VerbosityLevel > 0)  {
        G4cout << __func__ << " Activated step limit for volume "<<v->logical->GetName() << G4endl;
      }
    }
  }


  // helper function
  void Mu2eWorld::setStepLimitToAllSuchVolumes(const G4String& vn,
                                               G4UserLimits* const stepLimit,
                                               const G4LogicalVolumeStore* const lvs,
                                               int verbosityLevel) {
    int vtbcc = 0;
    for ( auto lvi=lvs->begin(); lvi!=lvs->end(); ++lvi) {
      if ((*lvi)->GetName() == vn) {
	(*lvi)->SetUserLimits( stepLimit );
        if(verbosityLevel > 0) {
          G4cout << __func__<< " Activated step limit for volume " << vn << G4endl;
        }
        ++vtbcc;
      }
    }
    if (vtbcc>1 && verbosityLevel > -1) {
      G4cout << __func__<< " WARNING: found " << vn << " " << vtbcc << " times" << G4endl;
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

    // G4LogicalVolume* hallAir        = _helper->locateVolInfo("HallAir").logical;

    G4LogicalVolume* ds1Vacuum      = _helper->locateVolInfo("DS1Vacuum").logical;
    G4LogicalVolume* ds2Vacuum      = _helper->locateVolInfo("DS2Vacuum").logical;
    G4LogicalVolume* ds3Vacuum      = _helper->locateVolInfo("DS3Vacuum").logical;
    G4LogicalVolume* dsCVacuum      = _helper->locateVolInfo("DSCryoVacuumRegion").logical;

    G4LogicalVolume* tracker        = _helper->locateVolInfo("TrackerMother").logical;
    G4LogicalVolume* trackerPanel   = _helper->locateVolInfo("StrawPanelEnvelope").logical;
    G4LogicalVolume* caloMother     = _helper->locateVolInfo("CalorimeterMother").logical;
    G4LogicalVolume* stoppingTarget = _helper->locateVolInfo("StoppingTargetMother").logical;
    G4LogicalVolume* productionTarget = _helper->locateVolInfo("ProductionTargetMother").logical;

    vector<G4LogicalVolume*> psVacua;
    psVacua.push_back( _helper->locateVolInfo("PSVacuum").logical );
    psVacua.push_back( _helper->locateVolInfo("psVacuumVesselVacuum").logical );

    vector<G4LogicalVolume*> tsVacua; // could use stepLimiterHelper; see below
    tsVacua.push_back( _helper->locateVolInfo("TS1Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS2Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS3Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS4Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS5Vacuum").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS1CryoInsVac").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS2CryoInsVac").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS3CryoInsVac").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS4CryoInsVac").logical );
    tsVacua.push_back( _helper->locateVolInfo("TS5CryoInsVac").logical );


    vector<G4LogicalVolume*> mbsLVS;
    mbsLVS.push_back( _helper->locateVolInfo("MBSMother").logical );

    // We may make separate G4UserLimits objects per logical volume but we choose not to.
    // At some it might be interesting to make several step limiters, each with different
    // limits.  For now that is not necessary.
    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();
    G4UserLimits* stepLimit = reg.add( G4UserLimits(bfieldMaxStep_) );
    if(_g4VerbosityLevel > 0) {
      G4cout << __func__ << " Using step limit = "<<bfieldMaxStep_/CLHEP::mm<<" mm"<<G4endl;
    }

    // Add the step limiters to the interesting volumes.
    // Keep them separated so that we can add different step limits should we decide to.

    // hallAir->SetUserLimits( stepLimit ); // not a vacuum per se; CPU costly

    // calorimeter does not use the nest functions; fixme
    // so we need to use the geant4's functions
    // in addition, some of its logical volumes are duplicated; fixme

    G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
    setStepLimitToAllSuchVolumes("caloBackPlateFEELog",    stepLimit,lvs,_g4VerbosityLevel);
    setStepLimitToAllSuchVolumes("caloCrystalFrameInLog",  stepLimit,lvs,_g4VerbosityLevel);
    setStepLimitToAllSuchVolumes("caloFEBLog",             stepLimit,lvs,_g4VerbosityLevel);
    setStepLimitToAllSuchVolumes("caloFEEBoxInLog",        stepLimit,lvs,_g4VerbosityLevel);
    setStepLimitToAllSuchVolumes("caloFrontPlateLog",      stepLimit,lvs,_g4VerbosityLevel);
    setStepLimitToAllSuchVolumes("caloFullBackPlateLog",   stepLimit,lvs,_g4VerbosityLevel);
    setStepLimitToAllSuchVolumes("caloHoleBackLog",        stepLimit,lvs,_g4VerbosityLevel);
    setStepLimitToAllSuchVolumes("calofullCrystalDiskLog", stepLimit,lvs,_g4VerbosityLevel);
    setStepLimitToAllSuchVolumes("ccrateBoxInLog",         stepLimit,lvs,_g4VerbosityLevel);
    setStepLimitToAllSuchVolumes("ccrateFullBoxLog",       stepLimit,lvs,_g4VerbosityLevel);

    ds1Vacuum->SetUserLimits( stepLimit );
    ds2Vacuum->SetUserLimits( stepLimit );
    ds3Vacuum->SetUserLimits( stepLimit );
    dsCVacuum->SetUserLimits( stepLimit );

    tracker->SetUserLimits( stepLimit );
    trackerPanel->SetUserLimits( stepLimit );

    caloMother->SetUserLimits( stepLimit );
    stoppingTarget->SetUserLimits( stepLimit );
    productionTarget->SetUserLimits( stepLimit );

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
    stepLimiterHelper("^TrackerPlaneEnvelope_.*$",                 stepLimit );
    stepLimiterHelper("^TrackerSupportServiceEnvelope_.*$",        stepLimit );
    stepLimiterHelper("^TrackerSupportServiceSectionEnvelope_.*$", stepLimit );

    // special case for straws to get the energy deposits right

    if (strawGasMaxStep_>0.0) {
      G4UserLimits* strawGasStepLimit = reg.add( G4UserLimits(strawGasMaxStep_) );
      if(_g4VerbosityLevel > 0) {
        G4cout << __func__<< " Using strawGas step limit = "
               <<strawGasMaxStep_/CLHEP::mm<<" mm"<<G4endl;
      }
      stepLimiterHelper("^TrackerStrawGas_.*$", strawGasStepLimit);
    }

    // and the calorimeter elements

    stepLimiterHelper("^caloDisk_.*$", stepLimit );
    stepLimiterHelper("^caloFEB_.*$", stepLimit );

    // calorimeter does not use Mu2e nest functions for those; see above; fixme
    // stepLimiterHelper("^caloFEBLog*$", stepLimit );

    // stepLimiterHelper("^ccrateFullBoxLog*$", stepLimit );
    // stepLimiterHelper("^ccrateBoxInLog*$", stepLimit );
    // stepLimiterHelper("^ccrateBoardLog*$", stepLimit );

    // stepLimiterHelper("^caloFrontPlateLog*$", stepLimit );
    // stepLimiterHelper("^caloCrystalFrameInLog*$", stepLimit );
    // stepLimiterHelper("^caloHoleBackLog*$", stepLimit );
    // stepLimiterHelper("^caloFEEBoxInLog*$", stepLimit );
    // stepLimiterHelper("^caloBackPlateFEELog*$", stepLimit );
    // stepLimiterHelper("^caloFullBackPlateLog*$", stepLimit );


    // An option to limit the step size in these non-vaccum volumes to
    // visually validate geometry of the filter channel
    //
    // ==== BEGIN COMMENT-OUT: to allow construction of new building dirt volumes w/o overlaps (knoepfel)
    //
    //     G4LogicalVolume* emfcMagnetAperture = _helper->locateVolInfo("ExtMonFNALfilterMagnetAperture").logical;
    //     if(emfcMagnetAperture) {
    //       const double maxStepLength = _config.getDouble("extMonFNAL.maxG4StepLength", 0)*CLHEP::millimeter;
    //       if(maxStepLength > 0) {
    //         G4cout<< __func__<<"Adding step limiter for ExtMonFNALFilterMagnet: maxStepLength = "<<maxStepLength<<G4endl;
    //         G4UserLimits* emfcStepLimit = reg.add(G4UserLimits(maxStepLength));
    //         emfcMagnetAperture->SetUserLimits(emfcStepLimit);
    //         _helper->locateVolInfo("ExtMonFNALfilterMagnetIron").logical->SetUserLimits(emfcStepLimit);
    //       }
    //     }
    //
    // ==== END COMMENT-OUT

    // Activate step limiter everywhere for spectial studies
    if(limitStepInAllVolumes_) {
      stepLimiterHelper("^.*$", stepLimit);
    }

  } // end Mu2eWorld::constructStepLimiters(){


  // Construct calorimeter if needed.
  VolumeInfo Mu2eWorld::constructCal(){

    // The calorimeter is built inside this volume.
    std::string theDS3("DS3Vacuum");
    if ( _config.getBool("inGaragePosition",false) ) theDS3 = "garageFakeDS3Vacuum";

    VolumeInfo const & detSolDownstreamVacInfo = _helper->locateVolInfo(theDS3);

    // Construct one of the calorimeters.
    VolumeInfo calorimeterInfo;
    if ( _config.getBool("hasDiskCalorimeter",false) ) {
      calorimeterInfo = constructDiskCalorimeter( detSolDownstreamVacInfo,_config );
    }

    return calorimeterInfo;

  }//Mu2eWorld::constructCal


  // A place holder for now.
  void Mu2eWorld::constructMagnetYoke(){
  }//Mu2eWorld::constructMagnetYoke


  // instantiateSensitiveDetectors
  void Mu2eWorld::instantiateSensitiveDetectors(){

      //art::ServiceHandle<GeometryService> geom;//I don't see this being used!
      G4SDManager* SDman = G4SDManager::GetSDMpointer();

      //get the LVStore singleton
      G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();

      // G4 takes ownership and will delete the detectors at the job end


/************************** Tracker **************************/
      //done
    if(sdHelper_->enabled(StepInstanceName::tracker)) {

        StrawSD* strawSD =
        new StrawSD( SensitiveDetectorName::TrackerGas(), _config );
        SDman->AddNewDetector(strawSD);

        //loop over all of the LV names and find the ones we need
        //set the SensitiveDetectors for these
        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
            G4String LVname = (*pos)->GetName();

            //from ConstructTrackerTDR and constructTrackerv3
            if (LVname.find("TrackerStrawGas_") != std::string::npos) {
              (*pos)->SetSensitiveDetector(strawSD);
            }
        }//for
     }//if tracker


/************************** TrackerDS **************************/
     //done
    if(sdHelper_->enabled(StepInstanceName::trackerDS)) {

        TrackerPlaneSupportSD* ttdsSD =
        new TrackerPlaneSupportSD( SensitiveDetectorName::TrackerPlaneSupport(), _config );
        SDman->AddNewDetector(ttdsSD);

        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
            G4String LVname = (*pos)->GetName();

            //from constructTrackerv3
            if (LVname.find("TrackerPlaneSupport_") != std::string::npos) {
              (*pos)->SetSensitiveDetector(ttdsSD);
            }

            //from constructTrackerv3Detailed
            if (LVname.find("TrackerSupportElecCu") != std::string::npos) {
              (*pos)->SetSensitiveDetector(ttdsSD);
            }
        }//for
    }//if trackerDS


/************************** VirtualDetector **************************/
    //done
    if(sdHelper_->enabled(StepInstanceName::virtualdetector)) {

        Mu2eSensitiveDetector* vdSD =
            new Mu2eSensitiveDetector( SensitiveDetectorName::VirtualDetector(), _config );
        SDman->AddNewDetector(vdSD);

        constructVirtualDetectorSDs(_config, vdSD);
    }//if virtualdetector


/************************** trackerSWires **************************/
    //done
    if (sdHelper_->enabled(StepInstanceName::trackerSWires)) {

        TrackerWireSD *ttwsSD =
            new TrackerWireSD(SensitiveDetectorName::TrackerSWires(), _config);
        SDman->AddNewDetector(ttwsSD);

        //trackerSWires need to be fixed if virtual detector is fixed
        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){

            G4String LVname = (*pos)->GetName();

            //from ConstructTrackerTDR
            if (LVname.find("TrackerStrawWire_") != std::string::npos) {
              (*pos)->SetSensitiveDetector(ttwsSD);
            }

            //from ConstructTrackerTDR
            if (LVname.find("TrackerWireCore_") != std::string::npos) {
              (*pos)->SetSensitiveDetector(ttwsSD);
            }

            //from ConstructTrackerTDR
            if (LVname.find("TrackerWirePlate_") != std::string::npos) {
              (*pos)->SetSensitiveDetector(ttwsSD);
            }
            }//for
    }//if trackerSWires


/************************** trackerWalls **************************/
    //done
    if (sdHelper_->enabled(StepInstanceName::trackerWalls)) {

        TrackerWireSD *ttwlSD =
            new TrackerWireSD( SensitiveDetectorName::TrackerWalls(), _config );
        SDman->AddNewDetector(ttwlSD);

        //trackerWalls need to be fixed if virtual detector is fixed
        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){

            G4String LVname = (*pos)->GetName();

            //from ConstructTrackerTDR and constructTrackerv3
            if (LVname.find("TrackerStrawWall_") != std::string::npos) {
              (*pos)->SetSensitiveDetector(ttwlSD);
            }

            //from ConstructTrackerTDR
            if (LVname.find("TrackerStrawWallOuterMetal_") != std::string::npos) {
              (*pos)->SetSensitiveDetector(ttwlSD);
            }

            //from ConstructTrackerTDR
            if (LVname.find("TrackerStrawWallInnerMetal1_") != std::string::npos) {
              (*pos)->SetSensitiveDetector(ttwlSD);
            }

            //from ConstructTrackerTDR
            if (LVname.find("TrackerStrawWallInnerMetal2_") != std::string::npos) {
              (*pos)->SetSensitiveDetector(ttwlSD);
            }
        }//for
    }//if trackerWalls


/************************** CALORIMETER **************************/
      // _geom is member data of Mu2eG4Universe, from which this inherits (GeoService const&)
      //done
    if ( const_cast<GeometryService&>(_geom).hasElement<Calorimeter>() ) {

      if(sdHelper_->enabled(StepInstanceName::calorimeter)) {
          CaloCrystalSD* ccSD =
            new CaloCrystalSD( SensitiveDetectorName::CaloCrystal(), _config );
          SDman->AddNewDetector(ccSD);

          for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
              G4String LVname = (*pos)->GetName();

              if (LVname.find("CrystalLog") != std::string::npos) {
                (*pos)->SetSensitiveDetector(ccSD);
              }
          }//for
      }//if calorimeter

      if(sdHelper_->enabled(StepInstanceName::calorimeterRO)) {
          CaloReadoutSD* crSD =
            new CaloReadoutSD( SensitiveDetectorName::CaloReadout(), _config );
          SDman->AddNewDetector(crSD);

          for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
              G4String LVname = (*pos)->GetName();

              if (LVname.find("CrystalROLog") != std::string::npos) {
                (*pos)->SetSensitiveDetector(crSD);
              }
          }//for
      }//if calorimeterRO

      if(sdHelper_->enabled(StepInstanceName::calorimeterROCard)) {
          CaloReadoutCardSD* crCardSD =
            new CaloReadoutCardSD( SensitiveDetectorName::CaloReadoutCard(), _config );
          SDman->AddNewDetector(crCardSD);

          for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
              G4String LVname = (*pos)->GetName();

              if (LVname.find("caloFEECardROLog") != std::string::npos) {
                (*pos)->SetSensitiveDetector(crCardSD);
              }
          }//for
      }//if calorimeterROCard

      if(sdHelper_->enabled(StepInstanceName::calorimeterCrate)) {
          CaloCrateSD* cCrateSD =
            new CaloCrateSD( SensitiveDetectorName::CaloCrate(), _config );
          SDman->AddNewDetector(cCrateSD);

          for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
              G4String LVname = (*pos)->GetName();

              if (LVname.find("ccrateActiveStripLog") != std::string::npos) {
                (*pos)->SetSensitiveDetector(cCrateSD);
              }
          }//for
      }//if calorimeterCrate

    }//if calorimeter


/************************** ExtMonFNALPixelSD **************************/
    if(true) { // this SD does not derive from Mu2eSensitiveDetector as it does not produce StepPointMCCollection
        GeomHandle<mu2e::ExtMonFNAL::ExtMon> extmon;
        //SDman->AddNewDetector(new ExtMonFNALPixelSD(_config, *extmon));

        ExtMonFNALPixelSD* emSD = new ExtMonFNALPixelSD(_config, *extmon);
        SDman->AddNewDetector(emSD);

        //from constructExtMonFNALModules
        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
            G4String LVname = (*pos)->GetName();

            if (LVname.find("EMFModule") != std::string::npos) {
              (*pos)->SetSensitiveDetector(emSD);
            }
        }//for
    }


/************************** StoppingTarget **************************/
    //done
    if(sdHelper_->enabled(StepInstanceName::stoppingtarget)) {

        Mu2eSensitiveDetector* stSD =
            new Mu2eSensitiveDetector( SensitiveDetectorName::StoppingTarget(), _config );
        SDman->AddNewDetector(stSD);

        //loop over all of the LV names to find ones we need
        //set the SensitiveDetectors for these
        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
            G4String LVname = (*pos)->GetName();

            if (LVname.find("Foil_") != std::string::npos) {
              (*pos)->SetSensitiveDetector(stSD);
            }

            if (LVname.find("FoilSupportStructure_") != std::string::npos) {
              (*pos)->SetSensitiveDetector(stSD);
            }
        }//for
    }//if stoppingtarget


/************************** CRV **************************/
    //done
    if(sdHelper_->enabled(StepInstanceName::CRV)) {

        Mu2eSensitiveDetector* sbSD =
            new Mu2eSensitiveDetector( SensitiveDetectorName::CRSScintillatorBar(), _config );
        SDman->AddNewDetector(sbSD);

        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
            G4String LVname = (*pos)->GetName();

            //we need names like CRV_R1, but not CRSCMB_CRV_R1 and not CRV_Support_TS
            if ((LVname.find("CRV_") != std::string::npos)
                && (LVname.find("CRV_Support") == std::string::npos)
                && (LVname.find("CRSCMB_CRV") == std::string::npos)) {
              (*pos)->SetSensitiveDetector(sbSD);
            }
        }//for
    }//if CRV


/************************** ProtonAbsorber **************************/
    if(sdHelper_->enabled(StepInstanceName::protonabsorber)) {
        Mu2eSensitiveDetector* paSD =
            new Mu2eSensitiveDetector( SensitiveDetectorName::ProtonAbsorber(),  _config );
        SDman->AddNewDetector(paSD);

        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
            G4String LVname = (*pos)->GetName();

            //from constructProtonAbsorber, will pick up protonabs{1,2,3,4}
            if (LVname.find("protonabs") != std::string::npos) {
              (*pos)->SetSensitiveDetector(paSD);
            }

            //from HelicalProtonAbsorber
            if (LVname.find("helical_pabs_log") != std::string::npos) {
              (*pos)->SetSensitiveDetector(paSD);
            }


        }//for
    }//if protonabsorber


/************************** PSVacuum **************************/
      //done
    if(sdHelper_->enabled(StepInstanceName::PSVacuum)) {

        Mu2eSensitiveDetector* psVacuumSD =
            new Mu2eSensitiveDetector( SensitiveDetectorName::PSVacuum(), _config );
        SDman->AddNewDetector(psVacuumSD);

        if( _config.getBool("PS.Vacuum.Sensitive", false) ) {
            psVacuumLogical_->SetSensitiveDetector(psVacuumSD);
        }
    }


/************************** STMDet **************************/
    //done
    if(sdHelper_->enabled(StepInstanceName::STMDet)) {
        Mu2eSensitiveDetector* STMDetSD =
        new Mu2eSensitiveDetector( SensitiveDetectorName::STMDet(), _config );
        SDman->AddNewDetector(STMDetSD);

        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){
            G4String LVname = (*pos)->GetName();

            if ( LVname == "stmDet1" || LVname == "stmDet2" ){
              (*pos)->SetSensitiveDetector(STMDetSD);
            }
        }//for
    }//if STMDet


/************************** panelEBKey **************************/
    if(sdHelper_->enabled(StepInstanceName::panelEBKey)) {
        Mu2eSensitiveDetector* EBKeySD =
        new Mu2eSensitiveDetector( SensitiveDetectorName::panelEBKey(), _config );
        SDman->AddNewDetector(EBKeySD);

        for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++) {
            G4String LVname = (*pos)->GetName();

            //we need the PanelEBKey, but NOT the PanelEBKeyShield
            if ((LVname.find("PanelEBKey") != std::string::npos) &&
                (LVname.find("PanelEBKeyShield") == std::string::npos)) {
              (*pos)->SetSensitiveDetector(EBKeySD);
            }
        }//for
    }//if panelEBKey


/************************** DSCableRun **************************/
      //if ( cableRunSensitive && sdHelper.enabled(StepInstanceName::DSCableRun) )
      if(sdHelper_->enabled(StepInstanceName::DSCableRun)) {
          Mu2eSensitiveDetector* cableRunSD =
          new Mu2eSensitiveDetector( SensitiveDetectorName::DSCableRun(), _config );
          SDman->AddNewDetector(cableRunSD);

          //NOTE: THIS 'if' test seems redundant to me, but I am just copying the format from constructDS.cc
          if ( _config.getBool("ds.CableRun.sensitive",false) ) {
              for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++) {
                  G4String LVname = (*pos)->GetName();

                  //will pick up names like CalCableRun, CalCableRunUpGap1, CalCableRunUpGap2, calCableRunFall
                  //but not CalCableRunLogInCalFeb, CalCableRunInCalFeb
                  if ( (LVname.find("alCableRun") != std::string::npos) &&
                       (LVname.find("CalCableRunLog") == std::string::npos) &&
                       (LVname.find("CalCableRunIn") == std::string::npos) )
                  {
                    (*pos)->SetSensitiveDetector(cableRunSD);
                  }

                  //will pick up names like TrkCableRun1, TrkCableRun2, TrkCableRunGap1, TrkCableRunGap1a, TrkCableRunGap2, TrkCableRunGap2a
                  //but not TrkCableRun1LogInCalFeb, TrkCableRun2LogInCalFeb, TrkCableRun1InCalFeb, TrkCableRun2InCalFeb
                  if ( (LVname.find("TrkCableRun") != std::string::npos) &&
                       (LVname.find("TrkCableRun1Log") == std::string::npos) &&
                       (LVname.find("TrkCableRun2Log") == std::string::npos) &&
                       (LVname.find("TrkCableRun1In") == std::string::npos) &&
                       (LVname.find("TrkCableRun2In") == std::string::npos) )
                  {
                    (*pos)->SetSensitiveDetector(cableRunSD);
                  }
              }//for
          }//if ds.CableRun.sensitive
      }//if DSCableRun

  }//instantiateSensitiveDetectors



} // end namespace mu2e
