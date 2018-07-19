#ifndef MCDataProducts_GenParticleCollections_hh
#define MCDataProducts_GenParticleCollections_hh

//
// Define a type for a collection of GenParticleCollections.
// Only used for use with the Tom LeCompte style G4MT
//
// Original author Rob Kutschke
//

#include <vector>

#include "MCDataProducts/inc/GenParticleCollection.hh"

namespace mu2e {
  typedef std::vector<mu2e::GenParticleCollection> GenParticleCollections;
}

#endif /* MCDataProducts_GenParticleCollections_hh */
