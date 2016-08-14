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
    explicit SimParticlesPrinter(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt) override;
  private:
    art::InputTag input_;
    SimParticleCollectionPrinter pr_;
  };

  //================================================================
  SimParticlesPrinter::SimParticlesPrinter(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , input_(pset.get<std::string>("inputCollection"))
    , pr_(pset.get<fhicl::ParameterSet>("printer", fhicl::ParameterSet()))
  {}

  //================================================================
  void SimParticlesPrinter::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<SimParticleCollection>(input_);
    pr_.print(std::cout, *ih);
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SimParticlesPrinter);
