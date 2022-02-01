//
//  This module creates the StrawGasStep objects used in downstream digi simulation, using the
//  G4 StepPointMCs
//
//  Original author: David Brown (LBNL), Krzysztof Genser 19 Aug. 2019
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

    for (const auto& step : *stepsHandle) {
      int tdc = 0;
      int adc = step.edep();
      STMDigi stm_digi(tdc, adc);
      outputSTMDigis->push_back(stm_digi);
    }

    event.put(std::move(outputSTMDigis));
  }
}

DEFINE_ART_MODULE(mu2e::MakeSTMDigis)
