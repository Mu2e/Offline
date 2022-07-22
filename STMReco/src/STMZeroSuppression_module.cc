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

// root
#include "TH1F.h"
#include "TF1.h"

#include "Offline/DataProducts/inc/STMTypes.hh"
#include "Offline/RecoDataProducts/inc/STMDigi.hh"
#include "Offline/RecoDataProducts/inc/STMHit.hh"
#include "Offline/STMReco/inc/ZPAlg.hh"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class STMZeroSuppression : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmDigisTag{ Name("stmDigisTag"), Comment("InputTag for STMDigiCollection")};
        fhicl::Atom<double> samplingFrequency{ Name("samplingFrequency"), Comment("Sampling Frequency of ADC [MHz]")};
        fhicl::Atom<double> tbefore{ Name("tbefore"), Comment("Store this time before the trigger [ns]")};
        fhicl::Atom<double> tafter{ Name("tafter"), Comment("Store this time after the trigger [ns]")};
        fhicl::Atom<double> threshold{ Name("threshold"), Comment("Threshold to define the trigger [ADC/ct]")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit STMZeroSuppression(const Parameters& conf);

    private:
    void beginJob() override;
      void produce(art::Event& e) override;

    art::InputTag _stmDigisTag;
    ZPAlg _zpAlg;
    double _ctPerNs;
  };

  STMZeroSuppression::STMZeroSuppression(const Parameters& config )  :
    art::EDProducer{config}
    ,_stmDigisTag(config().stmDigisTag())
    ,_zpAlg(config().samplingFrequency(), config().tbefore(), config().tafter(), config().threshold())
    ,_ctPerNs(1.0/(config().samplingFrequency()*1e-3))
  {
    consumes<STMDigiCollection>(_stmDigisTag);
    produces<STMDigiCollection>();
  }

  void STMZeroSuppression::beginJob() {
  }
    void STMZeroSuppression::produce(art::Event& event) {
    // create output
    unique_ptr<STMDigiCollection> outputSTMDigis(new STMDigiCollection);
    auto digisHandle = event.getValidHandle<STMDigiCollection>(_stmDigisTag);

    for (const auto& digi : *digisHandle) {
      std::cout << "Number of elements in unsuppressed file = " << digi.adcs().size() << std::endl;
      std::vector<size_t> starts;
      std::vector<size_t> ends;
      _zpAlg.ZeroSuppress(&digi.adcs()[0], digi.adcs().size(), starts, ends);
      for (size_t i = 0; i < starts.size(); ++i) {

        auto start = starts.at(i);
        auto end = ends.at(i);
        std::vector<int16_t> zp_adcs(digi.adcs().begin()+start, digi.adcs().begin()+end);
        STMDigi stm_digi(digi.trigNum(), STMTrigType(digi.trigType().mode(), digi.trigType().channel(), STMDataType::kZeroSuppressed), digi.trigTime()+_ctPerNs*start, 0, 0, 0, STMDigiFlag::kOK, zp_adcs);
        outputSTMDigis->push_back(stm_digi);
      }
    }

    event.put(std::move(outputSTMDigis));
  }
}

DEFINE_ART_MODULE(mu2e::STMZeroSuppression)
