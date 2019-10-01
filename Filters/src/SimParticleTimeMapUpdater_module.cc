// Update the map after compressing SimParticleCollection
//
// Andrei Gaponenko, 2014

#include <string>
#include <memory>

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"

namespace mu2e {

  class SimParticleTimeMapUpdater : public art::EDProducer {

  public:
    explicit SimParticleTimeMapUpdater(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& e) override;
  private:
    art::InputTag timeCollection_;
    art::InputTag remapping_;
  };

  //================================================================
  SimParticleTimeMapUpdater::SimParticleTimeMapUpdater(const fhicl::ParameterSet& pset)
    : art::EDProducer{pset}
    , timeCollection_(pset.get<std::string>("inputTimeCollection"))
    , remapping_(pset.get<std::string>("remapping"))
  {
    produces<SimParticleTimeMap>();
  }

  //================================================================
  void SimParticleTimeMapUpdater::produce(art::Event& event) {
    std::unique_ptr<SimParticleTimeMap> out(new SimParticleTimeMap);

    auto timeColl = event.getValidHandle<SimParticleTimeMap>(timeCollection_);
    auto remap = event.getValidHandle<SimParticleRemapping>(remapping_);

    for(const auto& entry : *timeColl) {
      const auto iter = remap->find(entry.first);
      if(iter != remap->end()) {
        (*out)[iter->second] = entry.second;
      }
    }

    event.put(std::move(out));
  }

  //================================================================
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::SimParticleTimeMapUpdater)
