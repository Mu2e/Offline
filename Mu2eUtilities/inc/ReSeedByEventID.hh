#ifndef Mu2eUtilities_ReSeedByEventID_hh
#define Mu2eUtilities_ReSeedByEventID_hh
//
// Reseed a random engine using the event id and a user supplied salt.
//

#include "Offline/SeedService/inc/SeedService.hh"

#include "canvas/Persistency/Provenance/EventID.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "CLHEP/Random/MixMaxRng.h"

namespace mu2e {

  class ReSeedByEventID{

  public:
    ReSeedByEventID( CLHEP::HepRandomEngine& engine, SeedService::seed_t salt);

    void reseed ( art::EventID const& id ) const;

  private:

    CLHEP::MixMaxRng *  _engine; // the engine to be reseeded
    SeedService::seed_t _salt;   // user supplied salt.

  };

} // namespace mu2e


#endif
