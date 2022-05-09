//
// Create the STMSteps from StepPointMCs
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include <utility>
// root
#include "TH1F.h"
#include "TTree.h"

#include "Offline/MCDataProducts/inc/STMStep.hh"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class MakeSTMSteps : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stepPointMCsTag{ Name("stepPointMCsTag"), Comment("InputTag for StepPointMCCollection")};
        fhicl::Atom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Level of verbosity to print")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit MakeSTMSteps(const Parameters& conf);

    private:
      void produce(art::Event& e) override;

    art::InputTag _stepPointMCsTag;
    int _verbosityLevel;
  };

  MakeSTMSteps::MakeSTMSteps(const Parameters& config )  :
    art::EDProducer{config}
    ,_stepPointMCsTag(config().stepPointMCsTag())
    ,_verbosityLevel(config().verbosityLevel())
  {
    consumes<StepPointMCCollection>(_stepPointMCsTag);
    produces<STMStepCollection>();
  }

  void MakeSTMSteps::produce(art::Event& event) {
    // create output
    unique_ptr<STMStepCollection> outputSTMSteps(new STMStepCollection);
    auto stepsHandle = event.getValidHandle<StepPointMCCollection>(_stepPointMCsTag);

    // For the time being just sum everything in the event FIXME
    if (stepsHandle->size() > 0) {
      double sum_edep = 0;
      for (const auto& step : *stepsHandle) {
        sum_edep += step.totalEDep();
      }
      double time = stepsHandle->at(0).time(); // just use the first step as a placeholder
      const auto& simPtr = stepsHandle->at(0).simParticle();
      STMStep stm_step(time, sum_edep, simPtr);
      outputSTMSteps->push_back(stm_step);
    }
    event.put(std::move(outputSTMSteps));
  }
}

DEFINE_ART_MODULE(mu2e::MakeSTMSteps)
