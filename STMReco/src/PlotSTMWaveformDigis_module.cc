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
        fhicl::Atom<std::string> xAxis{ Name("xAxis"), Comment("Choice of x-axis unit: \"sample_number\", \"waveform_time\", or \"event_time\"") };
        fhicl::Atom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Verbosity level")};
        fhicl::Atom<double> samplingFrequency{ Name("samplingFrequency"), Comment("Sampling Frequency of ADC [MHz]")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PlotSTMWaveformDigis(const Parameters& conf);

    private:
    void analyze(const art::Event& e) override;

    art::InputTag _stmWaveformDigisTag;
    bool _subtractPedestal;
    std::string _xAxis;
    int _verbosityLevel;
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;
    STMChannel _channel;
    double _ctPerNs;
  };

  PlotSTMWaveformDigis::PlotSTMWaveformDigis(const Parameters& config )  :
    art::EDAnalyzer{config},
    _stmWaveformDigisTag(config().stmWaveformDigisTag()),
    _subtractPedestal(config().subtractPedestal()),
    _xAxis(config().xAxis()),
    _verbosityLevel(config().verbosityLevel()),
    _ctPerNs((1.0/config().samplingFrequency())*1e3) // convert to ns
  {
    consumes<STMWaveformDigiCollection>(_stmWaveformDigisTag);
    _channel = STMUtils::getChannel(_stmWaveformDigisTag);
  }

  void PlotSTMWaveformDigis::analyze(const art::Event& event) {

    art::ServiceHandle<art::TFileService> tfs;
    auto waveformsHandle = event.getValidHandle<STMWaveformDigiCollection>(_stmWaveformDigisTag);

    std::stringstream histname, histtitle;
    int count = 0;
    STMEnergyCalib const& stmEnergyCalib = _stmEnergyCalib_h.get(event.id()); // get prodition
    const auto pedestal = stmEnergyCalib.pedestal(_channel);
    if (_verbosityLevel > 0) {
      std::cout << _channel.name() << " Pedestal = " << pedestal << std::endl;
    }

    double x_min = 0;
    double x_max = 0;
    int n_bins = 0;
    for (const auto& waveform : *waveformsHandle) {
      histname.str("");
      histname << "evt" << event.event() << "_waveform" << count;
      histtitle.str("");
      histtitle << "Event " << event.event() << " Waveform " << count << " (" << _channel.name() << ")";

      n_bins = waveform.adcs().size();
      if (_xAxis == "sample_number") {
        x_min = 0;
        x_max = n_bins;
      }
      else if (_xAxis == "waveform_time") {
        x_min = 0;
        x_max = n_bins*_ctPerNs;
      }
      else if (_xAxis == "event_time") {
        x_min = waveform.trigTimeOffset()*_ctPerNs;
        x_max = x_min+n_bins*_ctPerNs;
      }
      else {
        throw cet::exception("PlotSTMWaveformDigis") << "Invalid xAxis option: \"" << _xAxis << "\"" << std::endl;
      }

      TH1F* _hWaveform = tfs->make<TH1F>(histname.str().c_str(), histtitle.str().c_str(), n_bins,x_min,x_max);
      for (size_t i_adc = 0; i_adc < waveform.adcs().size(); ++i_adc) {
        const auto adc = waveform.adcs().at(i_adc);

        auto content = adc;
        if (_subtractPedestal) {
          content -= pedestal;
        }

        _hWaveform->SetBinContent(i_adc+1, content); // bins start numbering at 1
      }
      ++count;
    }
  }
}

DEFINE_ART_MODULE(mu2e::PlotSTMWaveformDigis)
