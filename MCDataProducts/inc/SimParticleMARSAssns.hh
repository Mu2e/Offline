// Andrei Gaponenko, 2012

#ifndef SimParticleMARSAssns_hh
#define SimParticleMARSAssns_hh

#include "canvas/Persistency/Common/Assns.h"

namespace mu2e {
  class SimParticle;
  class MARSInfo;
  typedef art::Assns<SimParticle, MARSInfo>  SimParticleMARSAssns;
}

#endif/*SimParticleMARSAssns_hh*/
