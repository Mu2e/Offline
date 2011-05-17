//
// Base class to allow generic access to all of the event generator classes.
// Also provides access to the random number engine associated with the 
// EventGenerator module.
//
// $Id: GeneratorBase.cc,v 1.2 2011/05/17 15:36:00 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:00 $
//
// Original author Rob Kutschke
// 

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"

namespace mu2e {

  art::RandomNumberGeneratorService::base_engine_t& 
  GeneratorBase::getEngine( const art::RandomNumberGeneratorService::label_t& engine_label ){

    static art::ServiceHandle<art::RandomNumberGeneratorService> rng;
    return rng->getEngine(engine_label);

  }

} // namespace mu2e
