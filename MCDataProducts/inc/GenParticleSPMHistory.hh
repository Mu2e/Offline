// Andrei Gaponenko, 2012

#ifndef GenParticleSPMHistory_hh
#define GenParticleSPMHistory_hh

#include "canvas/Persistency/Common/Assns.h"

namespace mu2e {
  class GenParticle;
  class StepPointMC;

  typedef art::Assns<GenParticle, StepPointMC>  GenParticleSPMHistory;
}

#endif/*GenParticleSPMHistory_hh*/
