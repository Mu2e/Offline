//
// Analyzer module to plot STM waveforms
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"

#include <utility>
#include <numeric>
// root
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"

#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/DataProducts/inc/STMChannel.hh"
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"

namespace mu2e {

  class PlotSTMWaveformDigis : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmWaveformDigisTag{ Name("stmWaveformDigisTag"), Comment("InputTag for STMWaveformDigiCollection")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PlotSTMWaveformDigis(const Parameters& conf);

    private:
    void analyze(const art::Event& e) override;

    art::InputTag _stmWaveformDigisTag;
    STMChannel _channel;
  };

  PlotSTMWaveformDigis::PlotSTMWaveformDigis(const Parameters& config )  :
    art::EDAnalyzer{config},
    _stmWaveformDigisTag(config().stmWaveformDigisTag())
  {
    consumes<STMWaveformDigiCollection>(_stmWaveformDigisTag);
    _channel = STMUtils::getChannel(_stmWaveformDigisTag);
  }

  void PlotSTMWaveformDigis::analyze(const art::Event& event) {

    art::ServiceHandle<art::TFileService> tfs;
    auto waveformsHandle = event.getValidHandle<STMWaveformDigiCollection>(_stmWaveformDigisTag);

    std::stringstream histname, histtitle;
    int count = 0;
    for (const auto& waveform : *waveformsHandle) {
      histname.str("");
      histname << "evt" << event.event() << "_waveform" << count;
      histtitle.str("");
      histtitle << "Event " << event.event() << " Waveform " << count << " (" << _channel.name() << ")";
      TH1F* _hWaveform = tfs->make<TH1F>(histname.str().c_str(), histtitle.str().c_str(), waveform.adcs().size(),0,waveform.adcs().size());
      for (size_t i_adc = 0; i_adc < waveform.adcs().size(); ++i_adc) {
        const auto adc = waveform.adcs().at(i_adc);
        _hWaveform->SetBinContent(i_adc+1, adc); // bins start numbering at 1
      }
      ++count;
    }
  }
}

DEFINE_ART_MODULE(mu2e::PlotSTMWaveformDigis)
