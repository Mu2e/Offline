//
// Primitive conditions data service.
// It does not yet do validty checking.
//
// $Id: ConditionsService.cc,v 1.7 2010/10/28 20:34:12 onoratog Exp $
// $Author: onoratog $
// $Date: 2010/10/28 20:34:12 $
//
// Original author Rob Kutschke
//

// C++ include files
#include <iostream>
#include <typeinfo>

// Framework include files
#include "DataFormats/Provenance/interface/ModuleDescription.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/Provenance/interface/RunID.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


// Mu2e include files
#include "ConditionsService/inc/ConditionsService.hh"

// Calibration entities.
// Would like to break the coupling to these.
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/PhysicsParams.hh"

using namespace std;

namespace mu2e {

  ConditionsService::ConditionsService(edm::ParameterSet const& iPS, 
                                       edm::ActivityRegistry&iRegistry) :
    _conditionsFile(iPS.getUntrackedParameter<std::string>("conditionsfile","conditions.txt")),
    _config(_conditionsFile),
    _entities(),
    _run_count()
  {
    iRegistry.watchPreBeginRun(this, &ConditionsService::preBeginRun);
  }
  

  ConditionsService::~ConditionsService(){
  }

  void 
  ConditionsService::preBeginRun(edm::RunID const& iID, edm::Timestamp const& iTime) {

    if(++_run_count > 1) {
      edm::LogWarning("CONDITIONS") << "This test version does not change geometry on run boundaries.";
      return;
    }
    
    edm::LogInfo log("CONDITIONS");
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
    addEntity( std::auto_ptr<ParticleDataTable>( new ParticleDataTable(_config)) );
    addEntity( std::auto_ptr<AcceleratorParams>( new AcceleratorParams(_config)) );
    addEntity( std::auto_ptr<DAQParams>        ( new DAQParams        (_config)) );
    addEntity( std::auto_ptr<PhysicsParams>    ( new PhysicsParams    (_config)) );
  }

  // Check that the configuration is self consistent.
  void ConditionsService::checkConsistency(){
  }

} // end namespace mu2e
