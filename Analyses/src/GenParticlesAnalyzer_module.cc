// Fills histograms for a GenParticle collection.
//
// Andrei Gaponenko, 2013

#include <iostream>
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

#include "Mu2eUtilities/inc/GeneratorSummaryHistograms.hh"

namespace mu2e {

  class GenParticlesAnalyzer : public art::EDAnalyzer {
  public:

    struct Config {
      fhicl::Atom<art::InputTag> inputs{
        fhicl::Name("inputs"),
          fhicl::Comment("The InputTag of a GenParticle collection to analyze. ")
          };
    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit GenParticlesAnalyzer(const Parameters& conf);

    void beginRun(const art::Run&) override;
    void analyze(const art::Event& e) override;

  private:
    art::InputTag inputs_;
    GeneratorSummaryHistograms genSummary_;
  };

  GenParticlesAnalyzer::GenParticlesAnalyzer(const Parameters& conf) :
      art::EDAnalyzer(conf),
      inputs_(conf().inputs()),
      genSummary_()
  {}

  // Can't book GeneratorSummaryHistograms in beginJob
  // because it uses the geometry service
  void GenParticlesAnalyzer::beginRun(const art::Run&) {
    genSummary_.book();
  }

  void GenParticlesAnalyzer::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<GenParticleCollection>(inputs_);
    genSummary_.fill(*ih);
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::GenParticlesAnalyzer);
