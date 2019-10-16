// Update StepPointMCs to point to a new SimParticle collection
//
// Andrei Gaponenko, 2016

#include <string>
#include <vector>
#include <memory>
#include <iostream>

#include "cetlib_except/exception.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "MCDataProducts/inc/SimParticleRemapping.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

namespace mu2e {

  class StepPointMCCollectionUpdater : public art::EDProducer {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> remapping {
        Name("remapping"),
          Comment("The SimParticleRemapping object to use.")
          };

      fhicl::Sequence<art::InputTag> inputs {
        Name("inputs"),
          Comment("A list of StepPointMCCollections to rewrite.")
          };

      fhicl::Atom<std::string> outInstanceName {
        Name("outInstanceName"),
          Comment(
                  "We create a single data product, so the normal Mu2e policy is\n"
                  "to not use instance name for it.  However when our output is\n"
                  "fed to FilterG4Out, the latter ignores module labels but keeps\n"
                  "instance names of its \"extraHitInputs\" collections.  If a job runs\n"
                  "several StepPointMCCollectionUpdater modules and their outputs\n"
                  "are passed to FilterG4Out, a way to preserve information is to\n"
                  "use instance names for StepPointMCCollectionUpdater outputs.\n"
                  )
          };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit StepPointMCCollectionUpdater(const Parameters& conf);

    void produce(art::Event& evt) override;
  private:
    typedef std::vector<art::InputTag> InputTags;
    art::InputTag remapping_;
    InputTags inputs_;
    std::string outInstanceName_;
  };

  //================================================================
  StepPointMCCollectionUpdater::StepPointMCCollectionUpdater(const Parameters& conf)
    : art::EDProducer{conf}
    , remapping_(conf().remapping())
    , inputs_(conf().inputs())
    , outInstanceName_(conf().outInstanceName())
  {
    produces<StepPointMCCollection>(outInstanceName_);
  }

  //================================================================
  void StepPointMCCollectionUpdater::produce(art::Event& event) {
    std::unique_ptr<StepPointMCCollection> out(new StepPointMCCollection());

    auto remap = event.getValidHandle<SimParticleRemapping>(remapping_);

    for(const auto& intag : inputs_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(intag);
      for(const auto& oldStep: *ih) {
        const auto oldPtr = oldStep.simParticle();
        const auto iter = remap->find(oldPtr);
        if(iter != remap->end()) {
          out->emplace_back(oldStep);
          out->back().simParticle() = iter->second;
        }
        else {
          throw cet::exception("BADCONFIG")
            <<"StepPointMCCollectionUpdater: SimParticleRemapping "<<remapping_
            <<" does not contain a particle used in StepPointMCCollection "<<intag
            <<"\n";
        }
      }
    }

    event.put(std::move(out), outInstanceName_);
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StepPointMCCollectionUpdater);
