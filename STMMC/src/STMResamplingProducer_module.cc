// Filters out the VD101 StepPointMCs ready for resampling
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
      fhicl::Atom<std::string> stepsTag{Name("VD101StepPointMCsTag"), Comment("Input tag of StepPointMCs associated with VD101")};
      fhicl::OptionalAtom<bool> verbose{Name("verbose"), Comment("Verbosity of output")};
    };

    using Parameters=art::EDProducer::Table<Config>;

    explicit STMResamplingProducer(const Parameters& pset);
    virtual void produce(art::Event& event) override;
    void beginJob() override;
    void endJob() override;

  private:
    art::InputTag _stepsTag;
    const uint16_t VirtualDetectorFilterID = 101; // Filter out all the StepPointMCs from VD101 for resampling
    bool _verbose;
    uint _keptStepPointMCCounter = 0;
  };
  // ===================================================
  STMResamplingProducer::STMResamplingProducer(const Parameters& config) :
    art::EDProducer{config},
    _stepsTag(config().stepsTag())
    {
      produces<StepPointMCCollection>();
      if (config().verbose.hasValue()) {_verbose = *std::move(config().verbose());}
      else {_verbose = false;}
    };
  // ===================================================
  void STMResamplingProducer::produce(art::Event& event)
  {
    // Define a handle to the virtualdetector
    art::Handle<StepPointMCCollection> _inputStepPointMCs;
    event.getByLabel(_stepsTag, _inputStepPointMCs);

    // Run check
    if (!(_inputStepPointMCs.isValid()))
      {
      if(_verbose){std::cout << "STMResamplingProducer : Handle is not valid, exiting." << std::endl;}
      exit(0);
      }
    if (_verbose){std::cout << "This event has " << _inputStepPointMCs->size() << " StepPointMCs with tag " << _stepsTag << std::endl;}

    // Define the StepPointMCCollection to be added to the event
    std::unique_ptr<StepPointMCCollection> _outputStepPointMCs(new StepPointMCCollection);

    // Check if the event has a hit in VirtualDetectorFilterID. If so add it to the collection
    for (const StepPointMC &step : *_inputStepPointMCs)
      {
      if ( step.volumeId() == VirtualDetectorFilterID){_outputStepPointMCs->emplace_back(step);};
      }

    _keptStepPointMCCounter += _outputStepPointMCs->size();
    event.put(std::move(_outputStepPointMCs));

    // return;
  };
  // ===================================================
  void STMResamplingProducer::beginJob()
  {
    if(_verbose){std::cout << "Looking for StepPointMCs with tag " << _stepsTag << std::endl;};
    // return;
  };
  // ===================================================
  void STMResamplingProducer::endJob()
  {
    if(_verbose){std::cout << "Number of kept StepPointMCs is " << _keptStepPointMCCounter << std::endl;};
    // return;
  };
}

DEFINE_ART_MODULE(mu2e::STMResamplingProducer)
