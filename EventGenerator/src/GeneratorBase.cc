//
// Base class to allow generic access to all of the event generator classes.
// Also provides access to the random number engine associated with the 
// EventGenerator module.
//
// $Id: GeneratorBase.cc,v 1.1 2010/08/18 06:31:19 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/08/18 06:31:19 $
//
// Original author Rob Kutschke
// 

// Framework includes
#include "FWCore/ServiceRegistry/interface/Service.h"

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"

namespace mu2e {

  edm::RandomNumberGeneratorService::base_engine_t& 
  GeneratorBase::getEngine( const edm::RandomNumberGeneratorService::label_t& engine_label ){

    static edm::Service<edm::RandomNumberGeneratorService> rng;
    return rng->getEngine(engine_label);

  }

} // namespace mu2e
