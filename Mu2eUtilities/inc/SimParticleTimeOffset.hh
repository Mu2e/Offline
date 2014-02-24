// A helper class to apply MC time offsets to account for e.g. proton
// pulse shape, or muon life time, to simulated particles.
//
// Andrei Gaponenko, 2014

#ifndef Mu2eUtilities_SimParticleTimeOffset_hh
#define Mu2eUtilities_SimParticleTimeOffset_hh

#include <vector>
#include <string>

#include "art/Utilities/InputTag.h"
#include "art/Persistency/Common/Ptr.h"

#include "MCDataProducts/inc/SimParticleTimeMap.hh"

namespace art { class Event; }
namespace fhicl { class ParameterSet; }

namespace mu2e {
  class StepPointMC;

  class SimParticleTimeOffset {
  public:
    SimParticleTimeOffset(const fhicl::ParameterSet& pset);

    void updateMap(const art::Event& evt);

    double totalTimeOffset(art::Ptr<SimParticle> p) const;
    double totalTimeOffset(const StepPointMC& s) const;
    double timeWithOffsetsApplied(const StepPointMC& s) const;

  private:
    std::vector<art::InputTag> inputs_;

    typedef std::vector<SimParticleTimeMap> Maps;
    mutable Maps offsets_;
  };
}

#endif/*Mu2eUtilities_SimParticleTimeOffset_hh*/
