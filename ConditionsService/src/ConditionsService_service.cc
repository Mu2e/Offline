//
// Primitive conditions data service.
// It does not yet do validty checking.
//
// $Id: ConditionsService_service.cc,v 1.7 2012/04/10 14:24:51 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/04/10 14:24:51 $
//
// Original author Rob Kutschke
//

// C++ include files
#include <iostream>
#include <typeinfo>

// Framework include files
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e include files
#include "ConditionsService/inc/ConditionsService.hh"

// Calibration entities.
// Would like to break the coupling to these.
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

using namespace std;

namespace mu2e {

  ConditionsService::ConditionsService(fhicl::ParameterSet const& pset,
                                       art::ActivityRegistry&iRegistry) :
    _conditionsFile(pset.get<std::string>("conditionsfile","conditions.txt")),
    _config(_conditionsFile,
            pset.get<bool>("allowReplacement",     true),
            pset.get<bool>("messageOnReplacement", true)
            ),
    _entities(),
    _run_count()
  {
    iRegistry.watchPreBeginRun(this, &ConditionsService::preBeginRun);
  }

  ConditionsService::~ConditionsService(){
  }

  // This template can be defined here because this is a private method which is only
  // used by the code below in the same file.
  template <typename ENTITY>
  void ConditionsService::addEntity(std::auto_ptr<ENTITY> d)
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

    mf::LogPrint log("CONDITIONS");
    log << "Conditions input file is: " << _conditionsFile << "\n";

    if ( _config.getBool("printConfig",false) ){
      log << "\n" << _config;
    }

    if ( _config.getBool("printConfigStats",false) ){
      // Work around absence of << operator for this print method.
      ostringstream os;
      _config.printStatistics(os);
      log << os.str();
    }

    checkConsistency();

    // Can we break the coupling to the entities?
    addEntity( std::auto_ptr<AcceleratorParams>  ( new AcceleratorParams  (_config)) );
    addEntity( std::auto_ptr<DAQParams>          ( new DAQParams          (_config)) );
    addEntity( std::auto_ptr<TrackerCalibrations>( new TrackerCalibrations(_config)) );
  }

  // Check that the configuration is self consistent.
  void ConditionsService::checkConsistency(){
  }

} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::ConditionsService);
