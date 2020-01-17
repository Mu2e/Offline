// Print out an EventIDSequence
//
// Andrei Gaponenko, 2020

#include <string>
#include <iostream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/IO/ProductMix/MixTypes.h"

#include "canvas/Utilities/InputTag.h"


namespace mu2e {

  //================================================================
  class EventIDSequencePrinter : public art::EDAnalyzer {
  public:


    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> input { Name("input"),
          Comment("The EventIDSequence to print") };

    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit EventIDSequencePrinter(const Parameters& conf);
    void analyze(const art::Event& evt) override;
  private:
    art::InputTag input_;
  };

  //================================================================
  EventIDSequencePrinter::EventIDSequencePrinter(const Parameters& pars)
    : art::EDAnalyzer(pars)
    , input_(pars().input())
  {}

  //================================================================
  void EventIDSequencePrinter::analyze(const art::Event& event) {
    auto ih = event.getValidHandle<art::EventIDSequence>(input_);
    std::cout<<"EventIDSequence "<<input_<<" : ";
    for(const auto& id : *ih) {
      std::cout<<" ( "<<id<<" ) ";
    }
    std::cout<<std::endl;
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EventIDSequencePrinter);
