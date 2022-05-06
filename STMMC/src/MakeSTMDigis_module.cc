//
// Create STMDigis from STMSteps
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
#include "TF1.h"

#include "Offline/MCDataProducts/inc/STMStep.hh"
#include "Offline/RecoDataProducts/inc/STMDigi.hh"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class MakeSTMDigis : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
       fhicl::Atom<art::InputTag> stmStepsTag{ Name("stmStepsTag"), Comment("InputTag for STMStepCollection")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit MakeSTMDigis(const Parameters& conf);

    private:
      void produce(art::Event& e) override;

    art::InputTag _stmStepsTag;
  };

  MakeSTMDigis::MakeSTMDigis(const Parameters& config )  :
    art::EDProducer{config},
    _stmStepsTag(config().stmStepsTag())
  {
    consumes<STMStepCollection>(_stmStepsTag);
    produces<STMDigiCollection>();
  }

    void MakeSTMDigis::produce(art::Event& event) {
    // create output
    unique_ptr<STMDigiCollection> outputSTMDigis(new STMDigiCollection);
    auto stepsHandle = event.getValidHandle<STMStepCollection>(_stmStepsTag);
    // Create Gaussian jitter
    TF1* res_fnc = new TF1("res_fnc", "TMath::Gaus(x, [0], [1])", -0.002, 0.002);
    res_fnc->SetParameters(0,0.001);

    for (const auto& step : *stepsHandle) {
      int tdc = 0;
      int adc = step.edep()*10000.0 + 10000.0*res_fnc->GetRandom();
      std::vector<int16_t> adcs;
      adcs.push_back(adc);
      STMDigi stm_digi(0,tdc,0,0,0,0, adcs);
      outputSTMDigis->push_back(stm_digi);
    }

    event.put(std::move(outputSTMDigis));
  }
}

DEFINE_ART_MODULE(mu2e::MakeSTMDigis)
