// Print out a StepPointMCCollection, based on StepPointsPrinter by 
// Andrei Gaponenko.
//
// In this version, can set a list of volumeIds to selectively print out.
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
  class SelectiveStepPtPrinter : public art::EDAnalyzer {
  public:
    explicit SelectiveStepPtPrinter(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt) override;
  private:
    std::string            input_;
    std::vector<int>       volumes_;
  };

  //================================================================
  SelectiveStepPtPrinter::SelectiveStepPtPrinter(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , input_(pset.get<std::string>("inputCollection"))
    , volumes_(pset.get<std::vector<int> >("volumes"))
  {}

  //================================================================
  void SelectiveStepPtPrinter::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<StepPointMCCollection>(input_);
    std::cout<<"Hits for "<<input_<<" in "<<event.id()<<": ("<<ih->size()<<")"<<std::endl;
    for(const auto& hit : *ih) {
      int id = hit.volumeId();
      CLHEP::Hep3Vector r = hit.position();
      for ( const auto& testId : volumes_ ) {
	//	std::cout << "id is " << id << ", while testId is " << testId << std::endl;
	if ( id == testId ) {
	  std::cout <<  "Volume:  " << id << ", particle: " << hit.simParticle()->pdgId() << ", position:  " << r << std::endl;
	}
      }
      //      std::cout<<"   "<<hit<<std::endl;
    }
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SelectiveStepPtPrinter);
