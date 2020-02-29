#include <exception>                            // for exception
#include <algorithm>                                   // for max
#include <map>                                         // for _Rb_tree_iterator
#include <memory>                                      // for allocator, uni...
#include <string>                                      // for string
#include <typeinfo>                                    // for type_info
#include <utility>                                     // for pair
#include <vector>                                      // for vector

#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"  // for SimParticleTim...
#include "cetlib_except/exception.h"                   // for exception, ope...
#include "fhiclcpp/ParameterSet.h"                     // for ParameterSet
#include "art/Framework/Principal/Event.h"             // for Event
#include "MCDataProducts/inc/StepPointMC.hh"           // for StepPointMC
#include "MCDataProducts/inc/StrawGasStep.hh"          // for StrawGasStep
#include "MCDataProducts/inc/SimParticle.hh"           // for SimParticle
#include "MCDataProducts/inc/SimParticleTimeMap.hh"    // for SimParticleTim...
#include "art/Framework/Principal/Handle.h"            // for ValidHandle
#include "canvas/Persistency/Common/Ptr.h"             // for Ptr, operator<<
#include "canvas/Utilities/InputTag.h"                 // for InputTag
#include "fhiclcpp/coding.h"                           // for ps_sequence_t
#include "fhiclcpp/exception.h"                        // for exception
#include "fhiclcpp/types/Sequence.h"                   // for Sequence


namespace mu2e {

  SimParticleTimeOffset::SimParticleTimeOffset(const Config& conf) {
    for(const auto& i: conf.inputs()) {
      inputs_.emplace_back(i);
    }
  }

  SimParticleTimeOffset::SimParticleTimeOffset(const fhicl::ParameterSet& pset) {
    // Do not add any offsets by default
    if(!pset.is_empty()) {
      typedef std::vector<std::string> VS;
      VS in(pset.get<VS>("inputs"));
      for(const auto& i : in) {
        inputs_.emplace_back(i);
      }
    }
  }

  SimParticleTimeOffset::SimParticleTimeOffset(const std::vector<art::InputTag>& tags) {
    for(const auto& i_tag : tags) {
      inputs_.emplace_back(i_tag);
    }
  }

  void SimParticleTimeOffset::updateMap(const art::Event& evt) {
    offsets_.clear();
    for(const auto& tag: inputs_) {
      auto m = evt.getValidHandle<SimParticleTimeMap>(tag);
      offsets_.emplace_back(*m);
    }
  }

  double SimParticleTimeOffset::totalTimeOffset(art::Ptr<SimParticle> p) const {

    if(offsets_.size() != inputs_.size()) {
      throw cet::exception("INVOCATION_ERROR")
        <<"SimParticleTimeOffset::totalTimeOffset():"
        <<" the number of loaded time maps "<<offsets_.size()
        <<" does not match the number of requested maps "<<inputs_.size()
        <<". Did you forget to call SimParticleTimeOffset::updateMap()?\n"
        ;
    }

    double dt = 0;

    // Look up the particle in all the maps, and add up the offsets
    for(auto& m : offsets_) {

      auto it = m.find(p);

      if(it == m.end()) { // no cached record for this particle

        const auto orig(p);

        // Navigate to the primary
        while(p->parent()) {
          p = p->parent();
        }

        it = m.find(p);
        if(it != m.end()) {
          // cache the result
          m[orig] = it->second;
        }
        else { // The ultimate parent must be in the map
          throw cet::exception("BADINPUTS")
            <<"SimParticleTimeOffset::totalTimeOffset(): the primary "<<p
            <<" is not in an input map\n";
        }
      } // caching

      dt += it->second;

    } // loop over offsets_ maps

    return dt;
  }

  double SimParticleTimeOffset::totalTimeOffset(const StepPointMC& s) const {
    return totalTimeOffset(s.simParticle());
  }

  double SimParticleTimeOffset::timeWithOffsetsApplied(const StepPointMC& s) const {
    return s.time() + totalTimeOffset(s);
  }

  double SimParticleTimeOffset::timeWithOffsetsApplied(const StrawGasStep& s) const {
    return s.time() + totalTimeOffset(s.simParticle());
  }

}
