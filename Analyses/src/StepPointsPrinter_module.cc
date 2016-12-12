// Print out a StepPointMCCollection.
//
// Andrei Gaponenko, 2013

#include <string>
#include <iostream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {

  //================================================================
  class StepPointsPrinter : public art::EDAnalyzer {
  public:
    explicit StepPointsPrinter(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt) override;
  private:
    std::string input_;
  };

  //================================================================
  StepPointsPrinter::StepPointsPrinter(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , input_(pset.get<std::string>("inputCollection"))
  {}

  //================================================================
  void StepPointsPrinter::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<StepPointMCCollection>(input_);
    std::cout<<"Hits for "<<input_<<" in "<<event.id()<<": ("<<ih->size()<<")"<<std::endl;
    for(const auto& hit : *ih) {
      std::cout<<"   "<<hit<<std::endl;
    }
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StepPointsPrinter);
