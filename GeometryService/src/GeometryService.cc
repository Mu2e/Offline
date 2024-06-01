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
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/DetectorSolenoidMaker.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/Mu2eHallMaker.hh"
#include "Offline/GeometryService/inc/TSdAMaker.hh"
#include "Offline/GeometryService/inc/TrackerMaker.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"
#include "Offline/GeometryService/inc/WorldG4Maker.hh"
#include "Offline/GeometryService/src/DetectorSystemMaker.hh"
#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Offline/ProductionTargetGeom/inc/ProductionTarget.hh"
#include "Offline/GeometryService/inc/ProductionTargetMaker.hh"
#include "Offline/ProductionSolenoidGeom/inc/ProductionSolenoid.hh"
#include "Offline/GeometryService/inc/ProductionSolenoidMaker.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "Offline/GeometryService/inc/PSEnclosureMaker.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "Offline/GeometryService/inc/PSVacuumMaker.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSShield.hh"
#include "Offline/GeometryService/inc/PSShieldMaker.hh"
#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "Offline/GeometryService/inc/ProtonBeamDumpMaker.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/GeometryService/inc/StoppingTargetMaker.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoidShielding.hh"
#include "Offline/GeometryService/inc/DetectorSolenoidShieldingMaker.hh"
#include "Offline/ExternalShieldingGeom/inc/ExtShieldUpstream.hh"
#include "Offline/GeometryService/inc/ExtShieldUpstreamMaker.hh"
#include "Offline/ExternalShieldingGeom/inc/ExtShieldDownstream.hh"
#include "Offline/GeometryService/inc/ExtShieldDownstreamMaker.hh"
#include "Offline/ExternalShieldingGeom/inc/Saddle.hh"
#include "Offline/GeometryService/inc/SaddleMaker.hh"
#include "Offline/ServicesGeom/inc/Pipe.hh"
#include "Offline/GeometryService/inc/PipeMaker.hh"
#include "Offline/ServicesGeom/inc/ElectronicRack.hh"
#include "Offline/GeometryService/inc/ElectronicRackMaker.hh"
#include "Offline/BeamlineGeom/inc/TSdA.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/DiskCalorimeterMaker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/BFieldGeom/inc/BFieldConfig.hh"
#include "Offline/GeometryService/inc/BFieldConfigMaker.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/BFieldManagerMaker.hh"
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/BeamlineGeom/inc/StraightSection.hh"
#include "Offline/GeometryService/inc/BeamlineMaker.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/GeometryService/inc/VirtualDetectorMaker.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/CosmicRayShieldMaker.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "Offline/GeometryService/inc/ExtMonFNALBuildingMaker.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "Offline/GeometryService/inc/ExtMonFNAL_Maker.hh"
#include "Offline/MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"
#include "Offline/GeometryService/inc/MECOStyleProtonAbsorberMaker.hh"
#include "Offline/MBSGeom/inc/MBS.hh"
#include "Offline/GeometryService/inc/MBSMaker.hh"
#include "Offline/STMGeom/inc/STM.hh"
#include "Offline/GeometryService/inc/STMMaker.hh"
#include "Offline/GeometryService/inc/Mu2eEnvelope.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMuonID.hh"
#include "Offline/GeometryService/inc/ExtMonFNALMuonIDMaker.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/GeometryService/inc/PTMMaker.hh"
#include "Offline/PTMGeom/inc/PTM.hh"

using namespace std;

namespace mu2e {

  GeometryService::GeometryService( const Parameters& pars,
                                    art::ActivityRegistry&iRegistry) :
    _inputfile(            pars().inputFile()),
    _bFieldFile(           pars().bFieldFile()),
    _allowReplacement(     pars().allowReplacement()),
    _messageOnReplacement( pars().messageOnReplacement()),
    _messageOnDefault(     pars().messageOnDefault()),
    _configStatsVerbosity( pars().configStatsVerbosity()),
    _printConfig(          pars().printConfig()),
    _printTopLevel(        pars().printConfigTopLevel()),
    _config(nullptr),
    _simulatedDetector(    pars.get_PSet().get<fhicl::ParameterSet>("simulatedDetector")),
    _standardMu2eDetector( _simulatedDetector.get<std::string>("tool_type") == "Mu2e"),
    _detectors()
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
      if( _run_count == 2 ){
        mf::LogWarning("GEOM") << "This test version does not change geometry on run boundaries.";
      }
      return;
    }

    _config = unique_ptr<SimpleConfig>(new SimpleConfig(_inputfile,
                                                      _allowReplacement,
                                                      _messageOnReplacement,
                                                      _messageOnDefault ));
    _config->printOpen(cout,"Geometry");

    _bfConfig = unique_ptr<SimpleConfig>(new SimpleConfig(_bFieldFile,
                                                          _allowReplacement,
                                                          _messageOnReplacement,
                                                          _messageOnDefault ));
    _bfConfig->printOpen(cout,"BField");


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
        //        std::cout << " adding Hayman in GeometryService" << std::endl;
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


    if (_config->getBool("hasPTM",false) ){
      std::unique_ptr<PTM> ptmon(PTMMaker::make(*_config));
      addDetector(std::move(ptmon));
    }

    if(_config->getBool("hasSTM",false)){
      STMMaker stm( *_config, beamline.solenoidOffset() );
      addDetector( stm.getSTMPtr() );
    }

    if(_config->getBool("hasVirtualDetector",false)){
      addDetector(VirtualDetectorMaker::make(*_config));
    }


    if(_bfConfig->getBool("hasBFieldManager",false)){
      std::unique_ptr<BFieldConfig> bfc( BFieldConfigMaker(*_bfConfig, beamline).getBFieldConfig() );
      BFieldManagerMaker bfmgr(*bfc);
      addDetector(std::move(bfc));
      addDetector(bfmgr.getBFieldManager());
    }


    if(_config->getBool("hasProtonAbsorber",false) && !_config->getBool("protonabsorber.isHelical", false) ){
      MECOStyleProtonAbsorberMaker mecopam( *_config, ds, target);
      addDetector( mecopam.getMECOStyleProtonAbsorberPtr() );
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
    if ( _configStatsVerbosity <= 0 ) return;
    ofstream gStats("geomStats.log");
    _config  ->printAllSummaries( gStats, _configStatsVerbosity, "" );
    ofstream bStats("bFieldStats.log");
    _bfConfig->printAllSummaries( bStats, _configStatsVerbosity, "" );
  }


} // end namespace mu2e
