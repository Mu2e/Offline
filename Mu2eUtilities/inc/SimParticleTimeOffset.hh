// A helper class to apply MC time offsets to account for e.g. proton
// pulse shape, or muon life time, to simulated particles.
//
// Andrei Gaponenko, 2014

#ifndef Mu2eUtilities_SimParticleTimeOffset_hh
#define Mu2eUtilities_SimParticleTimeOffset_hh

#include <vector>
#include <string>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h" // for legacy interface

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "MCDataProducts/inc/SimParticleTimeMap.hh"

namespace art { class Event; }
namespace fhicl { class ParameterSet; }

namespace mu2e {
  class StepPointMC;

  class SimParticleTimeOffset {
  public:

    struct Config {
      fhicl::Sequence<art::InputTag> inputs{
        fhicl::Name("inputs"),
          fhicl::Comment("List of SimParticleTimeMap collection tags to use."),
          std::vector<art::InputTag>()
          };
    };

    explicit SimParticleTimeOffset(const Config& conf);
    explicit SimParticleTimeOffset(const fhicl::ParameterSet& pset); // legacy
    explicit SimParticleTimeOffset(const std::vector<art::InputTag>& tags);

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
