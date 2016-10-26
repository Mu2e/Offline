// Andrei Gaponenko, 2012

#ifndef GenParticleMARSAssns_hh
#define GenParticleMARSAssns_hh

#include "canvas/Persistency/Common/Assns.h"

namespace mu2e {
  class GenParticle;
  class MARSInfo;
  typedef art::Assns<GenParticle, MARSInfo>  GenParticleMARSAssns;
}

#endif/*GenParticleMARSAssns_hh*/
