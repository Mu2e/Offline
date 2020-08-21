#ifndef EventGenerator_GeneratorBase_hh
#define EventGenerator_GeneratorBase_hh
//
// Base class to allow generic access to all of the event generator classes.
// Also provides access to the random number engine associated with the
// EventGenerator module.
//
//
// Original author Rob Kutschke
//

// Framework includes
#include "art/Framework/Core/detail/EngineCreator.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// Mu2e includes
#include "MCDataProducts/inc/GenParticleCollection.hh"

// CLHEP includes
#include "CLHEP/Random/RandomEngine.h"

namespace mu2e {

  class GeneratorBase{
  public:
    explicit GeneratorBase(bool const isFromG4bl = false)
      : _isFromG4bl{isFromG4bl}
    {}

    virtual ~GeneratorBase() = default;
    virtual void generate( GenParticleCollection&  ) = 0;

    bool isFromG4bl() { return _isFromG4bl; }

  protected:
    bool _isFromG4bl;
  };

} // end namespace mu2e,

#endif /* EventGenerator_GeneratorBase_hh */
