//
// Create STMHits from STMDigis
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

#include "Offline/RecoDataProducts/inc/STMDigi.hh"
#include "Offline/RecoDataProducts/inc/STMHit.hh"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class MakeSTMHits : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmDigisTag{ Name("stmDigisTag"), Comment("InputTag for STMDigiCollection")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit MakeSTMHits(const Parameters& conf);

    private:
      void produce(art::Event& e) override;

    art::InputTag _stmDigisTag;
  };

  MakeSTMHits::MakeSTMHits(const Parameters& config )  :
    art::EDProducer{config},
    _stmDigisTag(config().stmDigisTag())
  {
    consumes<STMDigiCollection>(_stmDigisTag);
    produces<STMHitCollection>();
  }

    void MakeSTMHits::produce(art::Event& event) {
    // create output
    unique_ptr<STMHitCollection> outputSTMHits(new STMHitCollection);
    auto digisHandle = event.getValidHandle<STMDigiCollection>(_stmDigisTag);

    for (const auto& digi : *digisHandle) {
      int tdc = digi.trigTime();
      int adc = digi.adcs().at(0); // pick the first ADC as a placeholder
      float time = tdc;
      float energy = adc/10000.0;
      STMHit stm_hit(time,energy);
      outputSTMHits->push_back(stm_hit);
    }

    event.put(std::move(outputSTMHits));
  }
}

DEFINE_ART_MODULE(mu2e::MakeSTMHits)
