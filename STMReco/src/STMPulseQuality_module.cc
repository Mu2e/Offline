//
// Create zero-suppressed STMDigis from unsuppressed STMDigis
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include <utility>
#include <algorithm>
#include <numeric>
// root
#include "TH1F.h"
#include "TF1.h"

#include "Offline/DataProducts/inc/STMTypes.hh"
#include "Offline/RecoDataProducts/inc/STMDigi.hh"
#include "Offline/STMReco/inc/PQAlg.hh"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class STMPulseQuality : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmDigisTag{ Name("stmDigisTag"), Comment("InputTag for STMDigiCollection")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit STMPulseQuality(const Parameters& conf);

    private:
    void beginJob() override;
    void produce(art::Event& e) override;

    art::InputTag _stmDigisTag;
    mu2e::PQAlg _pqAlg;
  };

  STMPulseQuality::STMPulseQuality(const Parameters& config )  :
    art::EDProducer{config}
    ,_stmDigisTag(config().stmDigisTag())
    ,_pqAlg()
  {
    consumes<STMDigiCollection>(_stmDigisTag);
    produces<STMDigiCollection>();
  }

  void STMPulseQuality::beginJob() {
  }
    void STMPulseQuality::produce(art::Event& event) {
    // create output
    unique_ptr<STMDigiCollection> outputSTMDigis(new STMDigiCollection);
    auto digisHandle = event.getValidHandle<STMDigiCollection>(_stmDigisTag);

    //    double threshold = -100;
    int count = 0;
    for (const auto& digi : *digisHandle) {
      //      double pedestal = std::reduce(digi.adcs().begin(), digi.adcs().begin()+100) / 100.0;
      _pqAlg.process_pulse(digi);
      size_t n_found_pulses = _pqAlg.getNFound();

      for (size_t i_pulse = 0; i_pulse < n_found_pulses; ++i_pulse) {
        std::vector<int16_t> pq_adcs;
        pq_adcs.push_back(_pqAlg.getEnergy(i_pulse));
        uint32_t extra(0);
        STMDigi stm_digi(STMTrigType(digi.trigType().mode(), digi.trigType().channel().id(), STMDataType::kPQ), digi.trigTime()*1e3, 0, extra, STMDigiFlag::kOK, pq_adcs);
        outputSTMDigis->push_back(stm_digi);
      }
      ++count;
    }
    std::cout << "AE: No. PQ Digis = " << outputSTMDigis->size() << std::endl;
    event.put(std::move(outputSTMDigis));
  }
}

DEFINE_ART_MODULE(mu2e::STMPulseQuality)
