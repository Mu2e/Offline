// This effectively forwards G4ParticleTable::FindParticle() so that
// other Mu2e packages can check whether a particle code will cause
// Geant to crash.   G4ParticleTable can not be used directly from
// other Mu2e packages; Mu2eG4 uses a special build environment to
// make it available.
//
// Andrei Gaponenko, 2012

#ifndef Mu2eG4_inc_checkParticlCodeForG4_hh
#define Mu2eG4_inc_checkParticlCodeForG4_hh

namespace mu2e {
  // returns true if the particle is known, false otherwise
  bool checkParticleCodeForG4(int pdgId);
}

#endif/*Mu2eG4_inc_checkParticlCodeForG4_hh*/
