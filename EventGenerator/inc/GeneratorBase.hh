#ifndef GENERATORBASE_HH
#define GENERATORBASE_HH
//
// Base class to allow generic access to all of the event generator classes.
// Also provides access to the random number engine associated with the 
// EventGenerator module.
//
// $Id: GeneratorBase.hh,v 1.3 2011/05/17 15:35:59 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:35:59 $
//
// Original author Rob Kutschke
// 

// Framework includes
#include "art/Framework/Core/RandomNumberGeneratorService.h"

// Mu2e includes
#include "ToyDP/inc/ToyGenParticleCollection.hh"

// CLHEP includes
#include "CLHEP/Random/RandomEngine.h"

namespace mu2e {

  class GeneratorBase{

  public:
    GeneratorBase( ){
    }
    virtual ~GeneratorBase(){
    }

    virtual void generate( ToyGenParticleCollection&  ) = 0;

  protected:

    // A helper function to access the random number engine associated with this module.
    static art::RandomNumberGeneratorService::base_engine_t& 
    getEngine( const art::RandomNumberGeneratorService::label_t& engine_label = "" );

};

} // end namespace mu2e,

#endif

