//
// Primitive conditions data service.
// It does not yet do validty checking.
//
// Original author Rob Kutschke
//

// C++ include files
#include <iostream>
#include <typeinfo>

// Framework include files
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e include files
#include "ConditionsService/inc/ConditionsService.hh"

// Calibration entities.
// Would like to break the coupling to these.
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/ExtMonFNALConditions.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
// sources of served objects
#include "TrackerConditions/inc/StrawElectronics.hh"
#include "TrackerConditions/inc/StrawPhysics.hh"

using namespace std;

namespace mu2e {

  ConditionsService::ConditionsService(fhicl::ParameterSet const& pset,
                                       art::ActivityRegistry&iRegistry) :
    _pset(pset),
    _conditionsFile(       pset.get<std::string>("conditionsfile","conditions.txt")),
    _allowReplacement(     pset.get<bool>       ("allowReplacement",      true)),
    _messageOnReplacement( pset.get<bool>       ("messageOnReplacement",  false)),
    _messageOnDefault(     pset.get<bool>       ("messageOnDefault",      false)),
    _configStatsVerbosity( pset.get<int>        ("configStatsVerbosity",  0)),
    _printConfig(          pset.get<bool>       ("printConfig",           false)),

    // FIXME: should be initialized in beginRun, not c'tor.
    _config(_conditionsFile, _allowReplacement, _messageOnReplacement, _messageOnDefault),
    _entities(),
    _run_count()
  {
    iRegistry.sPreBeginRun.watch(this, &ConditionsService::preBeginRun);
    iRegistry.sPostEndJob.watch (this, &ConditionsService::postEndJob   );
  }

  ConditionsService::~ConditionsService(){
  }

  // This template can be defined here because this is a private method which is only
  // used by the code below in the same file.
  template <typename ENTITY>
  void ConditionsService::addEntity(std::unique_ptr<ENTITY> d)
  {
    if(_entities.find(typeid(ENTITY).name())!=_entities.end())
      throw cet::exception("GEOM") << "failed to install conditions entity with type name "
                                   << typeid(ENTITY).name() << "\n";

    ConditionsEntityPtr ptr(d.release());
    _entities[typeid(ENTITY).name()] = ptr;
  }

  void
  ConditionsService::preBeginRun(art::Run const &) {

    if(++_run_count > 1) {
      mf::LogPrint("CONDITIONS") << "This test version does not change geometry on run boundaries.";
      return;
    }

    cout << "Conditions input file is: " << _conditionsFile << "\n";

    if ( _printConfig ){ _config.print(cout, "Conditions: "); }

    checkConsistency();

    // Can we break the coupling to the entities?
    std::unique_ptr<AcceleratorParams>  acctmp(new AcceleratorParams(_config));
    const AcceleratorParams& accp = *acctmp;
    addEntity( std::move(acctmp) );
    addEntity( std::move(std::unique_ptr<DAQParams>          ( new DAQParams          (_config))) );
    addEntity( std::move(std::unique_ptr<CalorimeterCalibrations>( new CalorimeterCalibrations(_config))) );
    addEntity( std::move(std::unique_ptr<TrackerCalibrations>( new TrackerCalibrations(_config))) );
    addEntity( std::move(std::unique_ptr<ExtMonFNALConditions>( new ExtMonFNALConditions(accp, _config))) );
    addEntity(std::move(std::unique_ptr<StrawElectronics>(new StrawElectronics(_pset.get<fhicl::ParameterSet>("StrawElectronics",fhicl::ParameterSet())))) );
    addEntity(std::move(std::unique_ptr<StrawPhysics>(new StrawPhysics(_pset.get<fhicl::ParameterSet>("StrawPhysics",fhicl::ParameterSet())))) );

  }

  // Check that the configuration is self consistent.
  void ConditionsService::checkConsistency(){
  }

  // Called after all modules have completed their end of job.
  void ConditionsService::postEndJob(){
    _config.printAllSummaries ( cout, _configStatsVerbosity, "Conditions: " );
  }

} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::ConditionsService);
