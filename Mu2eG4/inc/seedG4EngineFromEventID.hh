//
//  seedG4EngineFromEventID.hh
//
//  Seed the Geant4 random engine of the calling thread from an event ID.
//  The simulation of an event then depends only on its ID and the salt,
//  not on which events were simulated before it or on which thread.
//  Mu2eG4MT does this for every event; sequential Mu2eG4 does it when
//  its seedFromEventID parameter is true.
//

#ifndef Mu2eG4_seedG4EngineFromEventID_hh
#define Mu2eG4_seedG4EngineFromEventID_hh

#include "canvas/Persistency/Provenance/EventID.h"

#include <array>
#include <string>

namespace mu2e {

  // Returns the two seeds it set.
  std::array<long, 2> seedG4EngineFromEventID(const art::EventID& id, const std::string& salt);

}

#endif /* Mu2eG4_seedG4EngineFromEventID_hh */
