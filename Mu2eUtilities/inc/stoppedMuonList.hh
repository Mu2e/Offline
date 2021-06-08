// Extract a list of stopped muons from a SimParticleCollection.  The
// intended clients are stopped muon resampling generators, which
// operate on inputs pre-filtered in the previous simulation stage.
// So we are NOT doing things like finding muon stops in the stopping
// target here.  Instead we simply look for muons that have the
// muMinusCaptureAtRest stopping code, and return them all.
//
// Andrei Gaponenko, 2021

#ifndef Mu2eUtilities_inc_stoppedMuonList_hh
#define Mu2eUtilities_inc_stoppedMuonList_hh

#include <vector>

#include "art/Framework/Principal/Handle.h"

#include "MCDataProducts/inc/SimParticle.hh"

namespace mu2e {
  std::vector<art::Ptr<SimParticle> > stoppedMuonList(art::ValidHandle<SimParticleCollection> simh);
  std::vector<art::Ptr<SimParticle> > stoppedMuonList(art::Handle<SimParticleCollection> simh);
}

#endif/*Mu2eUtilities_inc_stoppedMuonList_hh*/
