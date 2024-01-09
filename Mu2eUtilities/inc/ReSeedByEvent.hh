#ifndef Mu2eUtilities_ReSeedByEvent_hh
#define Mu2eUtilities_ReSeedByEvent_hh
//
// Reseed a random engine using the event id.
//
//

#include "Offline/SeedService/inc/SeedService.hh"

#include "canvas/Persistency/Provenance/EventID.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "CLHEP/Random/MixMaxRng.h"

namespace mu2e {

  class ReSeedByEvent{

  public:
    ReSeedByEvent( CLHEP::HepRandomEngine& engine, SeedService::seed_t seed);

    void reseed ( art::EventID const& id ) const;

  private:

    CLHEP::MixMaxRng *  _engine;
    SeedService::seed_t _salt;

  };

} // namespace mu2e


#endif
