//
// Primitive conditions data service.
// It does not yet do validty checking.
//
// $Id: ConditionsService.cc,v 1.9 2011/05/17 15:35:59 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:35:59 $
//
// Original author Rob Kutschke
//

// C++ include files
#include <iostream>
#include <typeinfo>

// Framework include files
#include "art/Persistency/Provenance/ModuleDescription.h"
#include "art/Persistency/Provenance/EventID.h"
#include "art/Persistency/Provenance/Timestamp.h"
#include "art/Persistency/Provenance/SubRunID.h"
#include "art/Persistency/Provenance/RunID.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// Mu2e include files
#include "ConditionsService/inc/ConditionsService.hh"

// Calibration entities.
// Would like to break the coupling to these.
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/PhysicsParams.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

using namespace std;

namespace mu2e {

  ConditionsService::ConditionsService(fhicl::ParameterSet const& iPS, 
                                       art::ActivityRegistry&iRegistry) :
    _conditionsFile(iPS.get<std::string>("conditionsfile","conditions.txt")),
    _config(_conditionsFile),
    _entities(),
    _run_count()
  {
    iRegistry.watchPreBeginRun(this, &ConditionsService::preBeginRun);
  }
  

  ConditionsService::~ConditionsService(){
  }

  void 
  ConditionsService::preBeginRun(art::RunID const& iID, art::Timestamp const& iTime) {

    if(++_run_count > 1) {
      mf::LogWarning("CONDITIONS") << "This test version does not change geometry on run boundaries.";
      return;
    }
    
    mf::LogInfo log("CONDITIONS");
    log << "Conditions input file is: " << _conditionsFile << "\n";

    if ( _config.get<bool>("printConfig",false) ){
      log << "\n" << _config;
    }

    if ( _config.get<bool>("printConfigStats",false) ){
      // Work around absence of << operator for this print method.
      ostringstream os;
      _config.printStatistics(os);
      log << os.str();
    }

    checkConsistency();

    // Can we break the coupling to the entities?
    addEntity( std::auto_ptr<ParticleDataTable>( new ParticleDataTable(_config)) );
    addEntity( std::auto_ptr<AcceleratorParams>( new AcceleratorParams(_config)) );
    addEntity( std::auto_ptr<DAQParams>        ( new DAQParams        (_config)) );
    addEntity( std::auto_ptr<PhysicsParams>    ( new PhysicsParams    (_config)) );
    addEntity( std::auto_ptr<TrackerCalibrations>    ( new TrackerCalibrations    (_config)) );
  }

  // Check that the configuration is self consistent.
  void ConditionsService::checkConsistency(){
  }

} // end namespace mu2e
