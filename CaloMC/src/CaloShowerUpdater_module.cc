// Re-write an existing CaloShowerStepCollection updating SimParticle pointers.
// Need to run this if the original SimParticleCollection is filtered.
//
// A. Gaponenko, 2018

#include <iostream>
#include <string>

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace mu2e {

  class CaloShowerUpdater : public art::EDProducer {
  private:
    art::InputTag showerInput_;
    art::InputTag newSimParticles_;

  public:
    struct Config {
      fhicl::Atom<std::string> showerInput{ fhicl::Name("showerInput") };
      fhicl::Atom<std::string> newSimParticles{ fhicl::Name("newSimParticles") };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit CaloShowerUpdater(const Parameters& config) :
       art::EDProducer{config},
       showerInput_(config().showerInput()),
       newSimParticles_(config().newSimParticles())
    {
      produces<CaloShowerStepCollection>();
    }

    virtual void produce( art::Event& e) override;
  };

  //--------------------------------------------------------------------
  void  CaloShowerUpdater::produce(art::Event& event) {
    std::unique_ptr<CaloShowerStepCollection> out(new CaloShowerStepCollection);

    const auto& newParticles = event.getValidHandle<SimParticleCollection>(newSimParticles_);
    art::ProductID newParticlesPID = newParticles.id();
    const art::EDProductGetter *newParticlesGetter(event.productGetter(newParticlesPID));

    const auto& in = event.getValidHandle<CaloShowerStepCollection>(showerInput_);
    for(const auto& i: *in) {
      art::Ptr<SimParticle> oldPtr(i.simParticle());
      art::Ptr<SimParticle> newPtr(newParticlesPID, oldPtr->id().asUint(), newParticlesGetter);
      out->emplace_back(i);
      out->back().setSimParticle(newPtr);
    }

    event.put(std::move(out));
  }

}

using mu2e::CaloShowerUpdater;
DEFINE_ART_MODULE(CaloShowerUpdater);
