// Andrei Gaponenko, 2021

#include <string>
#include <iostream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"

namespace mu2e {

  //================================================================
  class PhysVolumeInfoPrinter : public art::EDAnalyzer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> subrunInput {Name("subrunInput"),
          Comment("Tag of a PhysicalVolumeInfoMultiCollection in a SubRun to print"),
          ""
          };

      fhicl::Atom<art::InputTag> eventInput {Name("eventInput"),
          Comment("Tag of a PhysicalVolumeInfoMultiCollection in an Event to print"),
          ""
          };

    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit PhysVolumeInfoPrinter(const Parameters& pars);

    void endSubRun(const art::SubRun& sr) override;
    void analyze(const art::Event& evt) override;

  private:
    art::InputTag subrunInput_;
    art::InputTag eventInput_;
  };

  //================================================================
  PhysVolumeInfoPrinter::PhysVolumeInfoPrinter(const Parameters& pars)
    : art::EDAnalyzer(pars)
    , subrunInput_(pars().subrunInput())
    , eventInput_(pars().eventInput())
  {}

  //================================================================
  void PhysVolumeInfoPrinter::analyze(const art::Event& event) {
    if(!eventInput_.empty()) {
      auto ih = event.getValidHandle<PhysicalVolumeInfoMultiCollection>(eventInput_);
      std::cout<<"PhysVolumeInfoPrinter: event "<<event.id()<<", tag "<<eventInput_<<", coll = "<<*ih<<std::endl;
    }
  }

  //================================================================
  void PhysVolumeInfoPrinter::endSubRun(const art::SubRun& sr) {
    if(!subrunInput_.empty()) {
      auto ih = sr.getValidHandle<PhysicalVolumeInfoMultiCollection>(subrunInput_);
      std::cout<<"PhysVolumeInfoPrinter: SubRun "<<sr.id()<<", tag "<<subrunInput_<<", coll = "<<*ih<<std::endl;
    }
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::PhysVolumeInfoPrinter);
