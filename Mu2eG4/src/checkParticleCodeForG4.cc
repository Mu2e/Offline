#include "Mu2eG4/inc/checkParticleCodeForG4.hh"

#include "Geant4/G4ParticleTable.hh"

namespace mu2e {
  bool checkParticleCodeForG4(int pdgId) {
    return G4ParticleTable::GetParticleTable()->FindParticle(pdgId);
  }
}
