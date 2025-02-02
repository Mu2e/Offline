// Used with STMResamplingFilter. Adds empty StepPointMC collections with the associated tag so the filter doesn't crash.
// Orginal author: Pawel Plesniak

// stdlib includes
#include <iostream>
#include <string>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/StepPointMC.hh"


namespace mu2e {
  class STMResamplingProducer : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stepPointMCsTag{Name("StepPointMCsTag"), Comment("Input tag of StepPointMCs")};
        fhicl::Atom<unsigned long> virtualDetectorID{Name("VirtualDetectorID"), Comment("Virtual detector ID to generate a separate StepPointMC collection from")};
      };
      using Parameters=art::EDProducer::Table<Config>;
      explicit STMResamplingProducer(const Parameters& pset);
      virtual void produce(art::Event& event) override;
      virtual void endJob() override;

    private:
      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      unsigned long virtualDetectorID = 0, includedStepPointMCs = 0;
  };

  STMResamplingProducer::STMResamplingProducer(const Parameters& conf) :
    art::EDProducer{conf},
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().stepPointMCsTag())),
    virtualDetectorID(conf().virtualDetectorID()) {
      produces<StepPointMCCollection>();
    };

  void STMResamplingProducer::produce(art::Event& event) {
    // Get the data product
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    if (StepPointMCs.size() == 0)
      throw cet::exception("DataError") << "Requested data not found";

    // Define the StepPointMCCollection to be added to the event
    std::unique_ptr<StepPointMCCollection> outputStepPointMCs(new StepPointMCCollection);

    // Collect valid data products
    for (const StepPointMC& step : StepPointMCs) {
      if (step.volumeId() == virtualDetectorID)
        outputStepPointMCs->emplace_back(step);
    };

    // Update counter
    includedStepPointMCs += outputStepPointMCs->size();
    // Add the new data products to the event
    event.put(std::move(outputStepPointMCs));
    return;
  };

  void STMResamplingProducer::endJob() {
    mf::LogInfo log("STMResamplingProducer summary");
    log << "No. kept StepPointMCs: " << includedStepPointMCs  << "\n";
    return;
  };
}; // namespace mu2e

DEFINE_ART_MODULE(mu2e::STMResamplingProducer)
