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

#include <utility>
#include <algorithm>

// root
#include "TH1F.h"
#include "TF1.h"

#include "Offline/RecoDataProducts/inc/STMWaveform.hh"
#include "Offline/STMReco/inc/ZPAlg.hh"

namespace mu2e {

  class STMPedestalSubtraction : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmWaveformsTag{ Name("stmWaveformsTag"), Comment("InputTag for STMWaveformCollection")};
        fhicl::Atom<int16_t> pedestal{Name("pedestal"), Comment("Pedestal value")}; // TODO: get from DB
        fhicl::Atom<int> verbosityLevel{Name("verbosityLevel"), Comment("Verbosity level")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit STMPedestalSubtraction(const Parameters& conf);

    private:
    void beginJob() override;
    void produce(art::Event& e) override;

    art::InputTag _stmWaveformsTag;
    int _verbosityLevel;
    int16_t _pedestal;
  };

  STMPedestalSubtraction::STMPedestalSubtraction(const Parameters& config )  :
    art::EDProducer{config}
    ,_stmWaveformsTag(config().stmWaveformsTag())
    ,_verbosityLevel(config().verbosityLevel())
    ,_pedestal(config().pedestal())
  {
    consumes<STMWaveformCollection>(_stmWaveformsTag);
    produces<STMWaveformCollection>();
  }

  void STMPedestalSubtraction::beginJob() {
  }

  void STMPedestalSubtraction::produce(art::Event& event) {
    // create output
    auto waveformsHandle = event.getValidHandle<STMWaveformCollection>(_stmWaveformsTag);
    unique_ptr<STMWaveformCollection> outputSTMWaveforms(new STMWaveformCollection());

    //    if (_verbosityLevel > 0) {
    //    }

    for (const auto& waveform : *waveformsHandle) {
      std::vector<int16_t> pedsub_adcs;
      pedsub_adcs.reserve(waveform.adcs().size());
      for (const auto& adc : waveform.adcs()) {
        // if we have truncated the pulse, then subtracting the pedestal will just roll us over
        if (adc - _pedestal > std::numeric_limits<int16_t>::min()) {
          pedsub_adcs.push_back(adc - _pedestal);
        }
        else {
          pedsub_adcs.push_back(std::numeric_limits<int16_t>::min());
        }
      }
      STMWaveform stm_waveform(waveform.trigTimeOffset(), pedsub_adcs);
      outputSTMWaveforms->push_back(stm_waveform);
    }

    if (_verbosityLevel > 0) {
      std::cout << outputSTMWaveforms->size() << " waveforms found" << std::endl;
    }
    event.put(std::move(outputSTMWaveforms));
  }
}

DEFINE_ART_MODULE(mu2e::STMPedestalSubtraction)
