// Used with STMResamplingFilter. Adds empty StepPointMC collections with the associated tag so the filter doesn't crash.
//
// Pawel Plesniak

// stdlib includes
#include <iostream>
#include <string>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "canvas/Utilities/InputTag.h"

// Offline includes
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

using namespace std;

namespace mu2e{
  class STMResamplingProducer : public art::EDProducer
  {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    struct Config
    {
      fhicl::Atom<std::string> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Input tag of StepPointMCs")};
      fhicl::Atom<uint> VirtualDetectorID{Name("VirtualDetectorID"), Comment("Virtual detector ID to filter events from")};
      fhicl::OptionalAtom<bool> verbose{Name("verbose"), Comment("Prints summary")};
    };
    using Parameters=art::EDProducer::Table<Config>;

    explicit STMResamplingProducer(const Parameters& pset);
    virtual void produce(art::Event& event) override;
    virtual void endJob() override;
  private:
    art::ProductToken<StepPointMCCollection> StepPointMCsToken;
    uint VD_ID = 0, eventsWithVD_IDStep = 0, eventsWithoutVD_IDStep = 0, keptStepPointMCsCount = 0;
    bool verbose = false, eventHasVD_IDSteps = false;
  };
  // ===================================================
  STMResamplingProducer::STMResamplingProducer(const Parameters& conf) :
    art::EDProducer{conf},
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    VD_ID(conf().VirtualDetectorID())
    {
      produces<StepPointMCCollection>();
      auto _verbose = conf().verbose();
      if(_verbose)verbose = *_verbose;
    };
  // ===================================================
  void STMResamplingProducer::produce(art::Event& event)
  {
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    eventHasVD_IDSteps = false;
    // Define the StepPointMCCollection to be added to the event
    std::unique_ptr<StepPointMCCollection> _outputStepPointMCs(new StepPointMCCollection);

    // Check if the event has a required data product. If so add it to the collection
    for (const StepPointMC& step : StepPointMCs)
      {
        if ( step.volumeId() == VD_ID){ // This can be a VirutalDetectorID.
          eventHasVD_IDSteps = true;
          keptStepPointMCsCount++;
          _outputStepPointMCs->emplace_back(step);
        };
      };

    if (eventHasVD_IDSteps == true)
      {
        eventsWithVD_IDStep++;
      }
    else
      {
        eventsWithoutVD_IDStep++;
      };

    event.put(std::move(_outputStepPointMCs));
    return;
  };
// ===================================================
  void STMResamplingProducer::endJob()
  {
    if (verbose)
      {
        std::cout << "Events with steps in VD"         << VD_ID << ": " << eventsWithVD_IDStep    << std::endl;
        std::cout << "Number of kept StepPointMCs in " << VD_ID << ": " << keptStepPointMCsCount  << std::endl;
        std::cout << "Events without steps in VD"      << VD_ID << ": " << eventsWithoutVD_IDStep << std::endl;
      };
    return;
  };
};

DEFINE_ART_MODULE(mu2e::STMResamplingProducer)
