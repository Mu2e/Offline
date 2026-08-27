#ifndef SeedService_PerEventReseed_hh
#define SeedService_PerEventReseed_hh
//
// Per-event deterministic reseeding of a CLHEP engine (Offline issue #849).
//
// Under Mu2eG4MT, events reach the legacy modules downstream of g4run in a
// scheduler-dependent order, so any module consuming a sequential random
// stream is irreproducible run-to-run.  G4 itself is immune because
// Mu2eG4WorkerRunManager reseeds per event from a hash of the EventID; this
// header applies the identical recipe to any seeded module.
//
// Gated on the environment variable MU2E_RESEED_PER_EVENT so that default
// behavior (and reproducibility of existing single-threaded productions) is
// completely unchanged unless explicitly enabled.
//
#include "CLHEP/Random/RandomEngine.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include <cstdlib>
#include <functional>
#include <string>

namespace mu2e {

  inline bool perEventReseedEnabled() {
    static const bool on = (std::getenv("MU2E_RESEED_PER_EVENT") != nullptr);
    return on;
  }

  // Same construction as Mu2eG4WorkerRunManager::initializeRun: hash a
  // string of run/subrun/event plus a per-module salt, mask to 31 bits.
  inline void perEventReseed(CLHEP::HepRandomEngine& engine,
                             art::EventID const& id,
                             std::string const& salt) {
    if (!perEventReseedEnabled()) return;
    std::string const msg = "r" + std::to_string(id.run())
                          + "s" + std::to_string(id.subRun())
                          + "e" + std::to_string(id.event()) + salt;
    std::hash<std::string> hf;
    engine.setSeed(static_cast<long>(hf(msg) & 0x7FFFFFFFUL), 0);
  }

} // namespace mu2e
#endif
