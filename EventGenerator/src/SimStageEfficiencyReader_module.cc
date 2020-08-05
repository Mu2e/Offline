// Reads the efficiency of a simulation stage
//
// Andrew Edmonds, 2020

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "cetlib_except/exception.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"

#include "MCDataProducts/inc/SimStageEfficiency.hh"

namespace mu2e {

  //================================================================
  class SimStageEfficiencyCalculator : public art::EDProducer {

  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> simStageEffTag{Name("simStageEffTag"), Comment("Input tag for mu2e::SimStageEfficiency")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit SimStageEfficiencyCalculator(const Parameters& conf);
    virtual void produce(art::Event& event) override;
    virtual void beginRun(art::Run& run) override;

  private:
    Config _conf;
  };

  //================================================================
  SimStageEfficiencyCalculator::SimStageEfficiencyCalculator(const Parameters& conf)
    : art::EDProducer{conf},
    _conf(conf())
  {
    produces<mu2e::SimStageEfficiency, art::InRun>();
  }

  //================================================================
  void SimStageEfficiencyCalculator::produce(art::Event& event) {
  }

  //================================================================
  void SimStageEfficiencyCalculator::beginRun(art::Run& run) {
    auto simStageEffH = run.getValidHandle<mu2e::SimStageEfficiency>(_conf.simStageEffTag());
    const auto* simStageEff = simStageEffH.product();
    std::cout << "Efficiency = " << simStageEff->numerator() << " / " << simStageEff->denominator() << " = " << simStageEff->efficiency() << std::endl;

    run.put(std::unique_ptr<SimStageEfficiency>(new SimStageEfficiency(*simStageEff)));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SimStageEfficiencyCalculator);
