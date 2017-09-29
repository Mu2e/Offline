//
// Maintain up to date geometry information and serve it to
// other services and to the modules.
//
// Original author Rob Kutschke
// Reco version adapted from original by David Norvil Brown (UofL)

// C++ include files
#include <iostream>
#include <utility>

// Framework include files
#include "canvas/Persistency/Provenance/ModuleDescription.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "canvas/Persistency/Provenance/SubRunID.h"
#include "canvas/Persistency/Provenance/RunID.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e include files
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "GeometryService/inc/RecoGeometryService.hh"
#include "GeometryService/inc/DetectorSolenoidMaker.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/Mu2eHallMaker.hh"
#include "GeometryService/inc/TTrackerMaker.hh"
#include "GeometryService/inc/WorldG4.hh"
//#include "GeometryService/inc/WorldG4Maker.hh"
#include "GeometryService/src/DetectorSystemMaker.hh"
//#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "GeometryService/inc/StoppingTargetMaker.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/DiskCalorimeterMaker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "GeometryService/inc/BFieldConfigMaker.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/BFieldManagerMaker.hh"
#include "STMGeom/inc/STM.hh"
#include "GeometryService/inc/STMMaker.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"

using namespace std;

namespace mu2e {

  RecoGeometryService::RecoGeometryService(fhicl::ParameterSet const& pset,
                                   art::ActivityRegistry&iRegistry) :
    _inputfile(            pset.get<std::string> ("inputFile",            "geom000.txt")),
    _allowReplacement(     pset.get<bool>        ("allowReplacement",     true)),
    _messageOnReplacement( pset.get<bool>        ("messageOnReplacement", false)),
    _messageOnDefault(     pset.get<bool>        ("messageOnDefault",     false)),
    _configStatsVerbosity( pset.get<int>         ("configStatsVerbosity", 0)),
    _printConfig(          pset.get<bool>        ("printConfig",          false)),
    _config(nullptr),
    standardMu2eDetector_( pset.get<bool>        ("standardMu2eDetector", false)),
    _detectors(),
    _run_count()
  {
    iRegistry.sPreBeginRun.watch(this, &RecoGeometryService::preBeginRun);
    iRegistry.sPostEndJob.watch (this, &RecoGeometryService::postEndJob );
  }

  RecoGeometryService::~RecoGeometryService(){
  }

  // This template can be defined here because this is a private method which is only
  // used by the code below in the same file.
  template <typename DET>
  void RecoGeometryService::addDetector(std::unique_ptr<DET> d)
  {
    if(_detectors.find(typeid(DET).name())!=_detectors.end()) {
      throw cet::exception("GEOM") << "failed to install detector with type name "
                                   << typeid(DET).name() << "\n";
    }

      DetectorPtr ptr(d.release());
      _detectors[typeid(DET).name()] = ptr;
  }

  template <typename DETALIAS, typename DET>
  void RecoGeometryService::addDetectorAliasToBaseClass(std::unique_ptr<DET> d)
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
  RecoGeometryService::preBeginRun(art::Run const &) {

    if(++_run_count > 1) {
      mf::LogWarning("GEOM") << "This test version does not change geometry on run boundaries.";
      return;
    }

    cout  << "Geometry input file is: " << _inputfile << "\n";

    _config = unique_ptr<SimpleConfig>(new SimpleConfig(_inputfile,
                                                      _allowReplacement,
                                                      _messageOnReplacement,
                                                      _messageOnDefault ));

    // Print final state of file after all substitutions.
    if ( _printConfig      ){ _config->print(cout, "Geom: ");       }

    // decide if this is standard Mu2e detector or something else ...

    if (!isStandardMu2eDetector() ||
        !_config->getBool("mu2e.standardDetector",false)) {
      cout  << "Non standard mu2e configuration, assuming it is intentional" << endl;
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


    std::unique_ptr<StoppingTarget> tmptgt(StoppingTargetMaker(detSys.getOrigin(), *_config).getTargetPtr());
    //    const StoppingTarget& target = *tmptgt.get();
    addDetector(std::move(tmptgt));

    if (_config->getBool("hasTTracker",false)){
      TTrackerMaker ttm( *_config );
      addDetector( ttm.getTTrackerPtr() );
    }


    if(_config->getBool("hasDiskCalorimeter",false)){
      DiskCalorimeterMaker calorm( *_config, _config->getDouble("mu2e.solenoidOffset") );
      addDetector( calorm.calorimeterPtr() );
      addDetectorAliasToBaseClass<Calorimeter>( calorm.calorimeterPtr() );  //add an alias to detector list
    }


    // if(_config->getBool("hasBFieldManager",false)){
    //   std::unique_ptr<BFieldConfig> bfc(BFieldConfigMaker(*_config, beamline).getBFieldConfig());
    //   BFieldManagerMaker bfmgr(*bfc);
    //   addDetector(std::move(bfc));
    //   addDetector(bfmgr.getBFieldManager());
    // }

  } // preBeginRun()

  // Check that the configuration is self consistent.
  void RecoGeometryService::checkConfig(){
    checkTrackerConfig();
  }

  // Check that the Tracker configuration is self consistent.
  void RecoGeometryService::checkTrackerConfig(){
    int ntrackers(0);
    string allTrackers;
    if ( _config->getBool("hasTTracker",false) ) {
      allTrackers += " TTracker";
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
  void RecoGeometryService::addWorldG4(const Mu2eHall& hall) {
    throw cet::exception("GEOM")
      << "Trying to call G4 Maker in reco-only job! \n";
  }

  // Called after all modules have completed their end of job.
  void   RecoGeometryService::postEndJob(){
    _config->printAllSummaries( cout, _configStatsVerbosity, "Geom: " );
  }



} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::RecoGeometryService);
