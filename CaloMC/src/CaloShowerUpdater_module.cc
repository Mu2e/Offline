// Re-write an existing CaloShowerStepCollection updating SimParticle pointers.
// Need to run this if the original SimParticleCollection is filtered.
//
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"

#include <iostream>
#include <string>

namespace mu2e {

  class CaloShowerUpdater : public art::EDProducer
  {
      public:
         struct Config
         {
             fhicl::Atom<std::string> showerInput     { fhicl::Name("showerInput")     };
             fhicl::Atom<std::string> newSimParticles { fhicl::Name("newSimParticles") };
         };

         explicit CaloShowerUpdater(const art::EDProducer::Table<Config>& config) :
            art::EDProducer{config},
            showerInput_(config().showerInput()),
            newSimParticles_(config().newSimParticles())
         {
            produces<CaloShowerStepCollection>();
         }

         virtual void produce( art::Event& e) override;


      private:
         art::InputTag showerInput_;
         art::InputTag newSimParticles_;
  };


  //--------------------------------------------------------------------
  void  CaloShowerUpdater::produce(art::Event& event)
  {
      std::unique_ptr<CaloShowerStepCollection> out(new CaloShowerStepCollection);

      const auto& newParticles       = event.getValidHandle<SimParticleCollection>(newSimParticles_);
      art::ProductID newParticlesPID = newParticles.id();
      const art::EDProductGetter *newParticlesGetter(event.productGetter(newParticlesPID));

      const auto& in = event.getValidHandle<CaloShowerStepCollection>(showerInput_);
      for (const auto& i: *in)
      {
         art::Ptr<SimParticle> oldPtr(i.simParticle());
         art::Ptr<SimParticle> newPtr(newParticlesPID, oldPtr->id().asUint(), newParticlesGetter);
         out->emplace_back(i);
         out->back().setSimParticle(newPtr);
      }

      event.put(std::move(out));
  }

}

using mu2e::CaloShowerUpdater;
DEFINE_ART_MODULE(CaloShowerUpdater)
