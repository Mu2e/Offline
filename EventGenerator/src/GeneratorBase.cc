//
// Base class to allow generic access to all of the event generator classes.
// Also provides access to the random number engine associated with the
// EventGenerator module.
//
// $Id: GeneratorBase.cc,v 1.3 2011/05/18 02:27:16 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:16 $
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
