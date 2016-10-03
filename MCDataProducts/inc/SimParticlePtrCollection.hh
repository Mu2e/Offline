#ifndef MCDataProducts_SimParticlePtrCollection_hh
#define MCDataProducts_SimParticlePtrCollection_hh

#include <vector>
#include "canvas/Persistency/Common/Ptr.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace mu2e {
  typedef std::vector<art::Ptr<mu2e::SimParticle> > SimParticlePtrCollection;
}

#endif /* MCDataProducts_SimParticlePtrCollection_hh */
