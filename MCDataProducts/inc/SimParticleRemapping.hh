// Once a SimParticleCollecton is compressed, art::Ptr<SimParticle>
// data members in other objects need to be updated.  This transient
// only class provides an {oldPtr => newPtr} mapping necessary for
// such updates.
//
// Andrei Gaponenko, 2014

#ifndef SimParticleRemapping_hh
#define SimParticleRemapping_hh

#include <map>
#include "canvas/Persistency/Common/Ptr.h"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace mu2e {
  typedef std::map<art::Ptr<SimParticle>, art::Ptr<SimParticle> >  SimParticleRemapping;
}

#endif/*SimParticleRemapping_hh*/
