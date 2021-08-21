// A data product to associate SimParticles to time offsets of their
// primary protons.  This is an ingredient to allow the proton pulse
// profile be applied after G4 simulations.
//
// Andrei Gaponenko, 2013

#ifndef SimParticleProtonPulseTimeAssns_hh
#define SimParticleProtonPulseTimeAssns_hh

#include <map>

#include "canvas/Persistency/Common/Ptr.h"

#include "Offline/MCDataProducts/inc/SimParticle.hh"

namespace mu2e {
  typedef std::map<art::Ptr<SimParticle>, double>  SimParticleTimeMap;
}

#endif/*SimParticleProtonPulseTimeAssns_hh*/
