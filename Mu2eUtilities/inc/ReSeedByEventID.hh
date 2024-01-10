#ifndef Mu2eUtilities_ReSeedByEventID_hh
#define Mu2eUtilities_ReSeedByEventID_hh
//
// Reseed a random engine using the event id and a user supplied salt.
//
// Does not work with all CLHEP::RandomEngine concrete types.
//   - see description in .cc file
//

#include "Offline/SeedService/inc/SeedService.hh"

#include "canvas/Persistency/Provenance/EventID.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"

namespace mu2e {

  class ReSeedByEventID{

  public:
    ReSeedByEventID( CLHEP::HepRandomEngine& engine, SeedService::seed_t salt, int verbosity=0 );

    void reseed ( art::EventID const& id ) const;

  private:

    CLHEP::HepRandomEngine& engine_;     // The engine to be reseeded
    SeedService::seed_t     salt_;       // User supplied salt.
    int                     verbosity_;  // Control amount of printout; 0=none.

  };

} // namespace mu2e

#endif
