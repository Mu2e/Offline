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
  };

  PlotSTMWaveformDigis::PlotSTMWaveformDigis(const Parameters& config )  :
    art::EDAnalyzer{config},
    _stmWaveformDigisTag(config().stmWaveformDigisTag())
  {
    consumes<STMWaveformDigiCollection>(_stmWaveformDigisTag);
  }

  void PlotSTMWaveformDigis::analyze(const art::Event& event) {

    art::ServiceHandle<art::TFileService> tfs;
    auto waveformsHandle = event.getValidHandle<STMWaveformDigiCollection>(_stmWaveformDigisTag);

    std::stringstream histname;
    int count = 0;
    for (const auto& waveform : *waveformsHandle) {
      histname.str("");
      histname << "evt" << event.event() << "_waveform" << count;
      TH1F* _hWaveform = tfs->make<TH1F>(histname.str().c_str(), "", waveform.adcs().size(),0,waveform.adcs().size());
      int i_bin = 1;
      for (const auto& adc : waveform.adcs()) {
        _hWaveform->SetBinContent(i_bin, adc);
        ++i_bin;
      }
      ++count;
    }
  }
}

DEFINE_ART_MODULE(mu2e::PlotSTMWaveformDigis)
