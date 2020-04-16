// Print out a SimParticleCollection.
//
// Andrei Gaponenko, 2013

#include <string>
#include <iostream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"

namespace mu2e {

  //================================================================
  class SimParticlesPrinter : public art::EDAnalyzer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> inputCollection {Name("inputCollection"),
          Comment("Collection to print")
          };

      fhicl::Table<SimParticleCollectionPrinter::Config> printer { Name("printer") };
    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit SimParticlesPrinter(const Parameters& pars);

    void analyze(const art::Event& evt) override;

  private:
    art::InputTag input_;
    SimParticleCollectionPrinter pr_;
  };

  //================================================================
  SimParticlesPrinter::SimParticlesPrinter(const Parameters& pars)
    : art::EDAnalyzer(pars)
    , input_(pars().inputCollection())
    , pr_(pars().printer())
  {}

  //================================================================
  void SimParticlesPrinter::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<SimParticleCollection>(input_);
    pr_.print(std::cout, *ih);
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SimParticlesPrinter);
