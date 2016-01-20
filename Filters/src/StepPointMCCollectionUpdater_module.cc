// Update StepPointMCs to point to a new SimParticle collection
//
// Andrei Gaponenko, 2016

#include <string>
#include <vector>
#include <memory>
#include <iostream>

#include "cetlib/exception.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "MCDataProducts/inc/SimParticleRemapping.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

namespace mu2e {

  class StepPointMCCollectionUpdater : public art::EDProducer {
  public:
    explicit StepPointMCCollectionUpdater(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt) override;
  private:
    typedef std::vector<art::InputTag> InputTags;
    art::InputTag remapping_;
    InputTags inputs_;

    // We create a single data product, so the normal Mu2e policy is
    // to not use instance name for it.  However when our output is
    // fed to FilterG4Out, the latter ignores module labels but keeps
    // instance names of "extraHitInputs" collections.  If a job runs
    // several StepPointMCCollectionUpdater modules and their outputs
    // are passed to FilterG4Out, a way to preserve information is to
    // use instance names for StepPointMCCollectionUpdater outputs.
    std::string outInstanceName_;
  };

  //================================================================
  StepPointMCCollectionUpdater::StepPointMCCollectionUpdater(const fhicl::ParameterSet& pset)
    : remapping_(pset.get<std::string>("remapping"))
    , inputs_(pset.get<std::vector<art::InputTag> >("inputs"))
    , outInstanceName_(pset.get<std::string>("outInstanceName"))
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
