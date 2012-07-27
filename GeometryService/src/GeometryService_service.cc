//
// Maintain up to date geometry information and serve it to
// other services and to the modules.
//
// $Id: GeometryService_service.cc,v 1.37 2012/07/27 19:42:31 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/27 19:42:31 $
//
// Original author Rob Kutschke
//

// C++ include files
#include <iostream>
#include <typeinfo>

// Framework include files
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Persistency/Provenance/ModuleDescription.h"
#include "art/Persistency/Provenance/EventID.h"
#include "art/Persistency/Provenance/Timestamp.h"
#include "art/Persistency/Provenance/SubRunID.h"
#include "art/Persistency/Provenance/RunID.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// Mu2e include files
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/src/DetectorSystemMaker.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/WorldG4Maker.hh"
#include "Mu2eBuildingGeom/inc/BuildingBasics.hh"
#include "Mu2eBuildingGeom/inc/BuildingBasicsMaker.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuildingMaker.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "ProductionTargetGeom/inc/ProductionTargetMaker.hh"
#include "ProductionSolenoidGeom/inc/ProductionSolenoid.hh"
#include "ProductionSolenoidGeom/inc/ProductionSolenoidMaker.hh"
#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "ProductionSolenoidGeom/inc/PSEnclosureMaker.hh"
#include "ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "ProductionSolenoidGeom/inc/PSVacuumMaker.hh"
#include "ProductionSolenoidGeom/inc/PSShield.hh"
#include "ProductionSolenoidGeom/inc/PSShieldMaker.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDumpMaker.hh"
#include "TargetGeom/inc/Target.hh"
#include "TargetGeom/inc/TargetMaker.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "LTrackerGeom/inc/LTrackerMaker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TTrackerGeom/inc/TTrackerMaker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "ITrackerGeom/inc/ITrackerMaker.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/CalorimeterMaker.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "BFieldGeom/inc/BFieldConfigMaker.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "BFieldGeom/inc/BFieldManagerMaker.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/BeamlineMaker.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/VirtualDetectorMaker.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShieldMaker.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuildingMaker.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL_Maker.hh"
#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCI.hh"
#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCIMaker.hh"
#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"
#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorberMaker.hh"
#include "MBSGeom/inc/MBS.hh"
#include "MBSGeom/inc/MBSMaker.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"

using namespace std;

namespace mu2e {

  GeometryService::GeometryService(fhicl::ParameterSet const& pset,
                                   art::ActivityRegistry&iRegistry) :
    _inputfile(            pset.get<std::string> ("inputFile",            "geom000.txt")),
    _allowReplacement(     pset.get<bool>        ("allowReplacement",     true)),
    _messageOnReplacement( pset.get<bool>        ("messageOnReplacement", false)),
    _messageOnDefault(     pset.get<bool>        ("messageOnDefault",     false)),
    _configStatsVerbosity( pset.get<int>         ("configStatsVerbosity", 0)),
    _printConfig(          pset.get<bool>        ("printConfig",          false)),
    _config(0),
    _detectors(),
    _run_count()
  {
    iRegistry.watchPreBeginRun(this, &GeometryService::preBeginRun);
    iRegistry.watchPostEndJob (this, &GeometryService::postEndJob );
  }

  GeometryService::~GeometryService(){
  }

  // This template can be defined here because this is a private method which is only
  // used by the code below in the same file.
  template <typename DET>
  void GeometryService::addDetector(std::auto_ptr<DET> d)
  {
    if(_detectors.find(typeid(DET).name())!=_detectors.end())
      throw cet::exception("GEOM") << "failed to install detector with type name "
                                   << typeid(DET).name() << "\n";

      DetectorPtr ptr(d.release());
      _detectors[typeid(DET).name()] = ptr;
  }

  void
  GeometryService::preBeginRun(art::Run const &) {

    if(++_run_count > 1) {
      mf::LogWarning("GEOM") << "This test version does not change geometry on run boundaries.";
      return;
    }

    cout  << "Geometry input file is: " << _inputfile << "\n";

    _config = auto_ptr<SimpleConfig>(new SimpleConfig(_inputfile,
                                                      _allowReplacement,
                                                      _messageOnReplacement,
                                                      _messageOnDefault ));

    // Print final state of file after all substitutions.
    if ( _printConfig      ){ _config->print(cout, "Geom: ");       }

    // Throw if the configuration is not self consistent.
    checkConfig();

    // This must be the first detector added since other makers may wish to use it.
    addDetector(DetectorSystemMaker::make( *_config));

    // Make a detector for every component present in the configuration.

    std::auto_ptr<Beamline> tmpBeamline(BeamlineMaker::make(*_config));
    const Beamline& beamline = *tmpBeamline.get();
    addDetector(tmpBeamline);

    std::auto_ptr<ProductionTarget> tmpProdTgt(ProductionTargetMaker::make(*_config, beamline.solenoidOffset()));
    const ProductionTarget& prodTarget = *tmpProdTgt.get();
    addDetector(tmpProdTgt);

    std::auto_ptr<ProductionSolenoid>
      tmpProductionSolenoid(ProductionSolenoidMaker(*_config, beamline.solenoidOffset()).getProductionSolenoidPtr());

    const ProductionSolenoid& ps = *tmpProductionSolenoid.get();
    addDetector(tmpProductionSolenoid);

    std::auto_ptr<PSEnclosure>
      tmpPSE(PSEnclosureMaker::make(*_config, ps.psEndRefPoint()));
    const PSEnclosure& pse = *tmpPSE.get();
    addDetector(tmpPSE);

    // The Z coordinate of the boundary between PS and TS vacua
    const double vacPS_TS_z = -beamline.getTS().torusRadius() - 2*beamline.getTS().getTS1().getHalfLength();
    addDetector(PSVacuumMaker::make(*_config, ps, pse, vacPS_TS_z));

    addDetector(PSShieldMaker::make(*_config, ps.psEndRefPoint(), prodTarget.position()));

    std::auto_ptr<BuildingBasics> tmpBasics(BuildingBasicsMaker::make(*_config));
    const BuildingBasics& buildingBasics = *tmpBasics.get();
    addDetector(tmpBasics);

    const double dumpFrontShieldingYmin = buildingBasics.detectorHallFloorTopY() - buildingBasics.detectorHallFloorThickness();
    const double dumpFrontShieldingYmax = dumpFrontShieldingYmin +
      buildingBasics.detectorHallInsideFullHeight() + buildingBasics.detectorHallCeilingThickness();

    std::auto_ptr<ProtonBeamDump> tmpDump(ProtonBeamDumpMaker::make(*_config, dumpFrontShieldingYmin, dumpFrontShieldingYmax));

    const ProtonBeamDump& dump = *tmpDump.get();
    addDetector(tmpDump);


    std::auto_ptr<Mu2eBuilding> tmpbld(Mu2eBuildingMaker::make(*_config, buildingBasics, dump));
    const Mu2eBuilding& building = *tmpbld.get();
    addDetector(tmpbld);

    if(_config->getBool("hasTarget",false)){
      TargetMaker targm( *_config );
      addDetector( targm.getTargetPtr() );
    }

    if(_config->getBool("hasLTracker",false)){
      LTrackerMaker ltm( *_config );
      addDetector( ltm.getLTrackerPtr() );
    } else if (_config->getBool("hasITracker",false)){
      ITrackerMaker itm( *_config );
      addDetector( itm.getITrackerPtr() );
    } else if (_config->getBool("hasTTracker",false)){
      TTrackerMaker ttm( *_config );
      addDetector( ttm.getTTrackerPtr() );
    }

    if(_config->getBool("hasCalorimeter",false)){
      CalorimeterMaker calorm( *_config, beamline.solenoidOffset() );
      addDetector( calorm.getCalorimeterPtr() );
    }

    if(_config->getBool("hasCosmicRayShield",false)){
      CosmicRayShieldMaker crs( *_config, beamline.solenoidOffset() );
      addDetector( crs.getCosmicRayShieldPtr() );
    }

    std::auto_ptr<ExtMonFNALBuilding> tmpemb(ExtMonFNALBuildingMaker::make(*_config, dump));
    const ExtMonFNALBuilding& emfb = *tmpemb.get();
    addDetector(tmpemb);
    if(_config->getBool("hasExtMonFNAL",false)){
      addDetector(ExtMonFNAL::ExtMonMaker::make(*_config, emfb));
    }

    if(_config->getBool("hasExtMonUCI",false)){
      ExtMonUCI::ExtMonMaker extmon( *_config );
      addDetector( extmon.getDetectorPtr() );
    }

    if(_config->getBool("hasMBS",false)){
      MBSMaker mbs( *_config, beamline.solenoidOffset() );
      addDetector( mbs.getMBSPtr() );
    }

    if(_config->getBool("hasVirtualDetector",false)){
      addDetector(VirtualDetectorMaker::make(*_config));
    }

    if(_config->getBool("hasBFieldManager",false)){
      std::auto_ptr<BFieldConfig> bfc(BFieldConfigMaker(*_config, beamline).getBFieldConfig());
      BFieldManagerMaker bfmgr(*bfc);
      addDetector(bfc);
      addDetector(bfmgr.getBFieldManager());
    }

    if(_config->getBool("hasProtonAbsorber",false) && !_config->getBool("protonabsorber.isHelical", false) ){
      MECOStyleProtonAbsorberMaker mecopam( *_config );
      addDetector( mecopam.getMECOStyleProtonAbsorberPtr() );
    }

    addDetector(std::auto_ptr<Mu2eEnvelope>(new Mu2eEnvelope(building, dump, emfb)));

  } // preBeginRun()

  // Check that the configuration is self consistent.
  void GeometryService::checkConfig(){
    checkTrackerConfig();
    checkCalorimeterConfig();
  }

  // Check that the Tracker configuration is self consistent.
  void GeometryService::checkTrackerConfig(){
    int ntrackers(0);
    string allTrackers;
    if ( _config->getBool("hasLTracker",false) ) {
      allTrackers += " LTracker";
      ++ntrackers;
    }
    if ( _config->getBool("hasTTracker",false) ) {
      allTrackers += " TTracker";
      ++ntrackers;
    }
    if ( _config->getBool("hasITracker",false) ) {
      allTrackers += " ITracker";
      ++ntrackers;
    }

    if ( ntrackers > 1 ){
      throw cet::exception("GEOM")
        << "This configuration has more than one tracker: "
        << allTrackers
        << "\n";
    }

  }

  // Check that the Tracker configuration is self consistent.
  void GeometryService::checkCalorimeterConfig(){
    int nCalorimeters(0);
    string allCalorimeters;
    if ( _config->getBool("hasCalorimeter",false) ) {
      allCalorimeters += " Calorimeter";
      ++nCalorimeters;
    }
    if ( _config->getBool("hasDiskCalorimeter",false) ) {
      allCalorimeters += " DiskCalorimeter";
      ++nCalorimeters;
    }

    if ( nCalorimeters > 1 ){
      throw cet::exception("GEOM")
        << "This configuration has more than one calorimeter: "
        << allCalorimeters
        << "\n";
    }

  }

  // WorldG4 could be added along with all the other detectors in preBeginRun().
  // However we don't want to make WorldG4 available in non-Geant jobs.
  // Therefore it is added by G4_module via this dedicated call.
  void GeometryService::addWorldG4() {
    addDetector(WorldG4Maker::make(*_config));
  }

  // Called after all modules have completed their end of job.
  void   GeometryService::postEndJob(){
    _config->printAllSummaries( cout, _configStatsVerbosity, "Geom: " );
  }



} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::GeometryService);
