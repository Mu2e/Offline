//
//  seedG4EngineFromEventID.cc
//
//  See the header.
//

#include "Offline/Mu2eG4/inc/seedG4EngineFromEventID.hh"

#include "Geant4/Randomize.hh"

#include <functional>

namespace mu2e {

  std::array<long, 2> seedG4EngineFromEventID(const art::EventID& id, const std::string& salt) {
    const std::string msg = "r" + std::to_string(id.run())
      + "s" + std::to_string(id.subRun())
      + "e" + std::to_string(id.event()) + salt;
    std::hash<std::string> hf;
    long seeds[3] = { static_cast<long>(hf(msg + "1") & 0xFFFFFFFF),
                      static_cast<long>(hf(msg + "2") & 0xFFFFFFFF),
                      0 };
    G4Random::setTheSeeds(seeds, -1);
    return { seeds[0], seeds[1] };
  }

}
