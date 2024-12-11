// Ed Callaghan
// Apply a functional prescale calculated from KalSeed by configurable rule
// September 2024

// stl
#include <memory>

// art
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/detail/EngineCreator.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"

// cetlib_except
#include "cetlib_except/exception.h"

// clhep
#include "CLHEP/Random/RandFlat.h"

// fhiclcpp
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/Blinding/inc/KalSeedPrescaleTool.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/SeedService/inc/SeedService.hh"

namespace mu2e{
  class KalSeedFunctionalPrescale: public art::EDProducer{
    public:
      struct Config{
        fhicl::Atom<art::InputTag> kalseed_tag{
          fhicl::Name("KalSeedCollection"),
          fhicl::Comment("art::InputTag of KalSeedCollection to prescale")
        };
        fhicl::DelegatedParameter tool_config{
          fhicl::Name("tool"),
          fhicl::Comment("Configuration for prescale calculation tool")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      KalSeedFunctionalPrescale(const Parameters&);

    protected:
      art::InputTag _kalseed_tag;
      std::unique_ptr<KalSeedPrescaleTool> _tool;
      art::RandomNumberGenerator::base_engine_t& _engine;
      std::unique_ptr<CLHEP::RandFlat> _uniform;

      bool Query(const KalSeed&);

    private:
      void produce(art::Event&);
  };

  KalSeedFunctionalPrescale::KalSeedFunctionalPrescale(const Parameters& config):
    art::EDProducer(config),
    _kalseed_tag(config().kalseed_tag()),
    _engine(createEngine(art::ServiceHandle<SeedService>()->getSeed())){
      // instantiate prescale calculator and rng
      auto tool_config = config().tool_config.get<fhicl::ParameterSet>();
      _tool = art::make_tool<KalSeedPrescaleTool>(tool_config);
      _uniform = std::make_unique<CLHEP::RandFlat>(_engine);
      // framework hook
      this->consumes<KalSeedCollection>(_kalseed_tag);
      this->produces<KalSeedCollection>();
    }

  bool KalSeedFunctionalPrescale::Query(const KalSeed& kalseed){
    // calculate probability to pass the track
    double acceptance_rate = _tool->AcceptanceRate(kalseed);
    double x = _uniform->fire();

    // pass the track with that probability
    auto rv = false;
    if (x < acceptance_rate){
      rv = true;
    }

    return rv;
  }

  void KalSeedFunctionalPrescale::produce(art::Event& event){
    auto handle = event.getValidHandle<KalSeedCollection>(_kalseed_tag);
    auto accepted = std::make_unique<KalSeedCollection>();
    for (const auto& kalseed: *handle){
      auto passed = this->Query(kalseed);
      if (passed){
        accepted->push_back(kalseed);
      }
    }

    event.put(std::move(accepted));
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::KalSeedFunctionalPrescale)
