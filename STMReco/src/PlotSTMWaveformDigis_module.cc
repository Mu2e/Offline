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
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"

namespace mu2e {

  class PlotSTMWaveformDigis : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmWaveformDigisTag{ Name("stmWaveformDigisTag"), Comment("InputTag for STMWaveformDigiCollection")};
        fhicl::Atom<bool> subtractPedestal{ Name("subtractPedestal"), Comment("True/False whether to subtract the pedestal before plotting")};
        fhicl::Atom<std::string> xAxis{ Name("xAxis"), Comment("Choice of x-axis unit: \"sample_number\", \"adcs_time\", or \"event_time\"") };
        fhicl::Atom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Verbosity level")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PlotSTMWaveformDigis(const Parameters& conf);

    private:
    void analyze(const art::Event& e) override;

    art::ProductToken<STMWaveformDigiCollection> _stmWaveformDigisToken;
    bool _subtractPedestal;
    std::string _xAxis;
    int _verbosityLevel;
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;
    STMChannel _channel;
  };

  PlotSTMWaveformDigis::PlotSTMWaveformDigis(const Parameters& config )  :
    art::EDAnalyzer{config},
    _stmWaveformDigisToken(consumes<STMWaveformDigiCollection>(config().stmWaveformDigisTag())),
    _subtractPedestal(config().subtractPedestal()),
    _xAxis(config().xAxis()),
    _verbosityLevel(config().verbosityLevel()),
    _channel(STMUtils::getChannel(config().stmWaveformDigisTag()))
  { }

  void PlotSTMWaveformDigis::analyze(const art::Event& event) {

    art::ServiceHandle<art::TFileService> tfs;
    auto waveformsHandle = event.getValidHandle(_stmWaveformDigisToken);

    std::stringstream histname, histtitle;
    int count = 0;
    STMEnergyCalib const& stmEnergyCalib = _stmEnergyCalib_h.get(event.id()); // get prodition
    const auto pedestal = stmEnergyCalib.pedestal(_channel);
    if (_verbosityLevel > 0) {
      std::cout << _channel.name() << " Pedestal = " << pedestal << std::endl;
    }

    const auto nsPerCt = stmEnergyCalib.nsPerCt(_channel);

    for (const auto& adcs : *waveformsHandle) {
      histname.str("");
      histname << "evt" << event.event() << "_adcs" << count;
      histtitle.str("");
      histtitle << "Event " << event.event() << " Waveform " << count << " (" << _channel.name() << ")";

      Binning binning = STMUtils::getBinning(adcs, _xAxis, nsPerCt);
      TH1F* hWaveform = tfs->make<TH1F>(histname.str().c_str(), histtitle.str().c_str(), binning.nbins(),binning.low(),binning.high());

      for (size_t i_adc = 0; i_adc < adcs.adcs().size(); ++i_adc) {
        const auto adc = adcs.adcs().at(i_adc);

        auto content = adc;
        if (_subtractPedestal) {
          content -= pedestal;
        }

        hWaveform->SetBinContent(i_adc+1, content); // bins start numbering at 1
      }
      ++count;
    }
  }
}

DEFINE_ART_MODULE(mu2e::PlotSTMWaveformDigis)
