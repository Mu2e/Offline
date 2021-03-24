//
// Maintain up to date geometry information and serve it to
// other services and to the modules.
//
// Original author Rob Kutschke
//

// C++ include files
#include <iostream>
#include <utility>

// Framework include files
#include "art/Persistency/Provenance/ModuleDescription.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "canvas/Persistency/Provenance/SubRunID.h"
#include "canvas/Persistency/Provenance/RunID.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e include files
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/DetectorSolenoidMaker.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/Mu2eHallMaker.hh"
#include "GeometryService/inc/TSdAMaker.hh"
#include "GeometryService/inc/TrackerMaker.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/WorldG4Maker.hh"
#include "GeometryService/src/DetectorSystemMaker.hh"
#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "GeometryService/inc/ProductionTargetMaker.hh"
#include "ProductionSolenoidGeom/inc/ProductionSolenoid.hh"
#include "GeometryService/inc/ProductionSolenoidMaker.hh"
#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "GeometryService/inc/PSEnclosureMaker.hh"
#include "ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "GeometryService/inc/PSVacuumMaker.hh"
#include "ProductionSolenoidGeom/inc/PSShield.hh"
#include "GeometryService/inc/PSShieldMaker.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "GeometryService/inc/ProtonBeamDumpMaker.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "GeometryService/inc/StoppingTargetMaker.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoidShielding.hh"
#include "GeometryService/inc/DetectorSolenoidShieldingMaker.hh"
#include "ExternalShieldingGeom/inc/ExtShieldUpstream.hh"
#include "GeometryService/inc/ExtShieldUpstreamMaker.hh"
#include "ExternalShieldingGeom/inc/ExtShieldDownstream.hh"
#include "GeometryService/inc/ExtShieldDownstreamMaker.hh"
#include "ExternalShieldingGeom/inc/Saddle.hh"
#include "GeometryService/inc/SaddleMaker.hh"
#include "ServicesGeom/inc/Pipe.hh"
#include "GeometryService/inc/PipeMaker.hh"
#include "ServicesGeom/inc/ElectronicRack.hh"
#include "GeometryService/inc/ElectronicRackMaker.hh"
#include "BeamlineGeom/inc/TSdA.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/DiskCalorimeterMaker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "GeometryService/inc/BFieldConfigMaker.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/BFieldManagerMaker.hh"
#include "BeamlineGeom/inc/Beamline.hh"
#include "BeamlineGeom/inc/StraightSection.hh"
#include "GeometryService/inc/BeamlineMaker.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/VirtualDetectorMaker.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "GeometryService/inc/CosmicRayShieldMaker.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "GeometryService/inc/ExtMonFNALBuildingMaker.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/ExtMonFNAL_Maker.hh"
#include "MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"
#include "GeometryService/inc/MECOStyleProtonAbsorberMaker.hh"
#include "MBSGeom/inc/MBS.hh"
#include "GeometryService/inc/MBSMaker.hh"
#include "STMGeom/inc/STM.hh"
#include "GeometryService/inc/STMMaker.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMuonID.hh"
#include "GeometryService/inc/ExtMonFNALMuonIDMaker.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

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
    _printTopLevel(        pset.get<bool>        ("printConfigTopLevel",  false)),
    _config(nullptr),
    _pset   (pset),
    standardMu2eDetector_( _pset.get<std::string>("simulatedDetector.tool_type") == "Mu2e"),
    _detectors(),
    _run_count()
  {
    iRegistry.sPreBeginRun.watch(this, &GeometryService::preBeginRun);
    iRegistry.sPostEndJob.watch (this, &GeometryService::postEndJob );
  }

  GeometryService::~GeometryService(){
  }

  // This template can be defined here because this is a private method which is only
  // used by the code below in the same file.
  template <typename DET>
  void GeometryService::addDetector(std::unique_ptr<DET> d)
  {
    if(_detectors.find(typeid(DET).name())!=_detectors.end()) {
      throw cet::exception("GEOM") << "failed to install detector with type name "
                                   << typeid(DET).name() << "\n";
    }

      DetectorPtr ptr(d.release());
      _detectors[typeid(DET).name()] = ptr;
  }

  template <typename DETALIAS, typename DET>
  void GeometryService::addDetectorAliasToBaseClass(std::unique_ptr<DET> d)
  {

        std::string OriginalName = typeid(DET).name();
        DetMap::iterator it(_detectors.find(OriginalName));

        if(it==_detectors.end())
          throw cet::exception("GEOM")
            << "Can not alias an inexistant detector, detector " << OriginalName << "\n";

        std::string detectorName= typeid(DETALIAS).name() ;
        _detectors[detectorName] = it->second;
  }

  void
  GeometryService::preBeginRun(art::Run const &) {

    if(++_run_count > 1) {
      mf::LogWarning("GEOM") << "This test version does not change geometry on run boundaries.";
      return;
    }

    _config = unique_ptr<SimpleConfig>(new SimpleConfig(_inputfile,
                                                      _allowReplacement,
                                                      _messageOnReplacement,
                                                      _messageOnDefault ));

    _config->printOpen(cout,"Geometry");

    if(_printTopLevel) {
      //print the top level geometry file contents
      //the top level often contains a single named config file or a list of specific version files
      ConfigFileLookupPolicy configFile;
      std::string file = configFile(_inputfile);
      std::ifstream in(file.c_str());
      if ( !in ) {
	// No conf file for this test.
	throw cet::exception("Geom")
	  << "GeometryService: Cannot open input file: "
	  << file
	  << endl;
      }
      std::cout << "GeometryService: printing top level geometry file:\n";
      std::string line;
      while ( in ){
	std::getline(in,line);
	if ( !in ){
	  break;
	}

	std::cout << line.c_str() << std::endl;
      }
      std::cout << "GeometryService: finished printing top level geometry file.\n";
    }
    // Print final state of file after all substitutions.
    if ( _printConfig      ){ _config->print(cout, "Geom: ");       }

    // 2019-03-24 P.M. : *not needed* decide if this is standard Mu2e detector or something else ...

    if (!isStandardMu2eDetector() ||
	!_config->getBool("mu2e.standardDetector",true)) {
      cout  << "Non-standard mu2e configuration, assuming it is intentional" << endl;
      return;
    }

    // Initialize geometry options
    _g4GeomOptions = unique_ptr<G4GeometryOptions>( new G4GeometryOptions( *_config ) );

    // Throw if the configuration is not self consistent.
    checkConfig();

    // This must be the first detector added since other makers may wish to use it.
    std::unique_ptr<DetectorSystem> tmpDetSys(DetectorSystemMaker::make(*_config));
    const DetectorSystem& detSys = *tmpDetSys.get();
    addDetector(std::move(tmpDetSys));

    // Make a detector for every component present in the configuration.

    std::unique_ptr<Beamline> tmpBeamline(BeamlineMaker::make(*_config));
    const Beamline& beamline = *tmpBeamline.get();
    addDetector(std::move(tmpBeamline));

    std::unique_ptr<ProductionTarget> tmpProdTgt(ProductionTargetMaker::make(*_config, beamline.solenoidOffset()));
    const ProductionTarget& prodTarget = *tmpProdTgt.get();
    addDetector(std::move(tmpProdTgt));

    std::unique_ptr<ProductionSolenoid>
      tmpProductionSolenoid(ProductionSolenoidMaker(*_config, beamline.solenoidOffset()).getProductionSolenoidPtr());

    const ProductionSolenoid& ps = *tmpProductionSolenoid.get();
    addDetector(std::move(tmpProductionSolenoid));

    std::unique_ptr<PSEnclosure>
      tmpPSE(PSEnclosureMaker::make(*_config, ps.psEndRefPoint()));
    const PSEnclosure& pse = *tmpPSE.get();
    addDetector(std::move(tmpPSE));

    // The Z coordinate of the boundary between PS and TS vacua
    StraightSection const * ts1vac = beamline.getTS().getTSVacuum<StraightSection>( TransportSolenoid::TSRegion::TS1 );
    const double vacPS_TS_z = ts1vac->getGlobal().z() - ts1vac->getHalfLength();

    addDetector(PSVacuumMaker::make(*_config, ps, pse, vacPS_TS_z));

    //addDetector(PSShieldMaker::make(*_config, ps.psEndRefPoint(), prodTarget.position()));

   if (_config->getString("targetPS_model") == "MDC2018"){ 
     //      std::cout << "adding Tier1 in GeometryService" << std::endl;
      addDetector(PSShieldMaker::make(*_config, ps.psEndRefPoint(), prodTarget.position()));
	} else 
      if (_config->getString("targetPS_model") == "Hayman_v_2_0"){ 
	//	std::cout << " adding Hayman in GeometryService" << std::endl;
	addDetector(PSShieldMaker::make(*_config, ps.psEndRefPoint(), prodTarget.haymanProdTargetPosition()));
	  } else 
	{throw cet::exception("GEOM") << " " << __func__ << " illegal production target version specified in GeometryService_service = " << _config->getString("targetPS_model")  << std::endl;}






    // Construct building solids
    std::unique_ptr<Mu2eHall> tmphall(Mu2eHallMaker::makeBuilding(*_g4GeomOptions,*_config));
    const Mu2eHall& hall = *tmphall.get();

    // Determine Mu2e envelope from building solids
    std::unique_ptr<Mu2eEnvelope> mu2eEnv (new Mu2eEnvelope(hall,*_config));

    // Make dirt based on Mu2e envelope
    Mu2eHallMaker::makeDirt( *tmphall.get(), *_g4GeomOptions, *_config, *mu2eEnv.get() );
    Mu2eHallMaker::makeTrapDirt( *tmphall.get(), *_g4GeomOptions, *_config, *mu2eEnv.get() );

    addDetector(std::move( tmphall ) );
    addDetector(std::move( mu2eEnv ) );

    std::unique_ptr<ProtonBeamDump> tmpDump(ProtonBeamDumpMaker::make(*_config, hall));
    const ProtonBeamDump& dump = *tmpDump.get();
    addDetector(std::move(tmpDump));

    // beamline info used to position DS
    std::unique_ptr<DetectorSolenoid> tmpDS( DetectorSolenoidMaker::make( *_config, beamline ) );
    const DetectorSolenoid& ds = *tmpDS.get();
    addDetector(std::move(tmpDS));

    // DS info used to position DS downstream shielding
    addDetector( DetectorSolenoidShieldingMaker::make( *_config, ds ) );

    std::unique_ptr<StoppingTarget> tmptgt(StoppingTargetMaker(detSys.getOrigin(), *_config).getTargetPtr());
    const StoppingTarget& target = *tmptgt.get();
    addDetector(std::move(tmptgt));

    if (_config->getBool("hasTracker",false)){
      TrackerMaker ttm( *_config );
      addDetector( ttm.getTrackerPtr() );
    }

    if(_config->getBool("hasMBS",false)){
      MBSMaker mbs( *_config, beamline.solenoidOffset() );
      addDetector( mbs.getMBSPtr() );
    }


    if(_config->getBool("hasDiskCalorimeter",false)){
      DiskCalorimeterMaker calorm( *_config, beamline.solenoidOffset() );
      addDetector( calorm.calorimeterPtr() );
      addDetectorAliasToBaseClass<Calorimeter>( calorm.calorimeterPtr() );  //add an alias to detector list
    }

    if(_config->getBool("hasCosmicRayShield",false)){
      CosmicRayShieldMaker crs( *_config, beamline.solenoidOffset() );
      addDetector( crs.getCosmicRayShieldPtr() );
    }

    if(_config->getBool("hasTSdA",false)){
      addDetector( TSdAMaker::make(*_config,ds) );
    }

    if(_config->getBool("hasExternalShielding",false)) {
      addDetector( ExtShieldUpstreamMaker::make(*_config)  );
      addDetector( ExtShieldDownstreamMaker::make(*_config));
      addDetector( SaddleMaker::make(*_config));
      addDetector( PipeMaker::make(*_config));
      addDetector( ElectronicRackMaker::make(*_config));
    }



    std::unique_ptr<ExtMonFNALBuilding> tmpemb(ExtMonFNALBuildingMaker::make(*_config, hall, dump));
    const ExtMonFNALBuilding& emfb = *tmpemb.get();
    addDetector(std::move(tmpemb));
    if(_config->getBool("hasExtMonFNAL",false)){
      addDetector(ExtMonFNAL::ExtMonMaker::make(*_config, emfb));
      addDetector(ExtMonFNALMuonIDMaker::make(*_config));
    }
    


    if(_config->getBool("hasVirtualDetector",false)){
      addDetector(VirtualDetectorMaker::make(*_config));
    }
      
    
    if(_config->getBool("hasBFieldManager",false)){
      std::unique_ptr<BFieldConfig> bfc( BFieldConfigMaker(*_config, beamline).getBFieldConfig() );
      BFieldManagerMaker bfmgr(*bfc);
      addDetector(std::move(bfc));
      addDetector(bfmgr.getBFieldManager());
    }
 

    if(_config->getBool("hasProtonAbsorber",false) && !_config->getBool("protonabsorber.isHelical", false) ){
      MECOStyleProtonAbsorberMaker mecopam( *_config, ds, target);
      addDetector( mecopam.getMECOStyleProtonAbsorberPtr() );
    }

    if(_config->getBool("hasSTM",false)){
      STMMaker stm( *_config, beamline.solenoidOffset() );
      addDetector( stm.getSTMPtr() );
    }

  } // preBeginRun()

  // Check that the configuration is self consistent.
  void GeometryService::checkConfig(){
    checkTrackerConfig();
  }

  // Check that the Tracker configuration is self consistent.
  void GeometryService::checkTrackerConfig(){
    int ntrackers(0);
    string allTrackers;
    if ( _config->getBool("hasTracker",false) ) {
      allTrackers += " Tracker";
      ++ntrackers;
    }
    if ( ntrackers > 1 ){
      throw cet::exception("GEOM")
        << "This configuration has more than one tracker: "
        << allTrackers
        << "\n";
    }

  }


  // WorldG4 could be added along with all the other detectors in preBeginRun().
  // However we don't want to make WorldG4 available in non-Geant jobs.
  // Therefore it is added by G4_module via this dedicated call.
  void GeometryService::addWorldG4(const Mu2eHall& hall) {
    addDetector(WorldG4Maker::make(hall,*_config));
  }

  // Called after all modules have completed their end of job.
  void   GeometryService::postEndJob(){
    _config->printAllSummaries( cout, _configStatsVerbosity, "Geom: " );
  }



} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::GeometryService);
