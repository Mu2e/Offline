//
// Base class to allow generic access to all of the event generator classes.
// Also provides access to the random number engine associated with the
// EventGenerator module.
//
// $Id: GeneratorBase.cc,v 1.4 2011/05/18 16:21:55 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/18 16:21:55 $
//
// Original author Rob Kutschke
//

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"

// Mu2e includes
#include "EventGenerator/inc/GeneratorBase.hh"

namespace mu2e {

  art::RandomNumberGenerator::base_engine_t&
  GeneratorBase::getEngine( const art::RandomNumberGenerator::label_t& engine_label ){

    static art::ServiceHandle<art::RandomNumberGenerator> rng;
    return rng->getEngine(engine_label);

  }

} // namespace mu2e
