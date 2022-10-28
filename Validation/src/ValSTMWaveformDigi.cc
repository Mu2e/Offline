#include "Offline/Validation/inc/ValSTMWaveformDigi.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include <algorithm>

int mu2e::ValSTMWaveformDigi::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _hNwf = tfs.make<TH1D>("Nwf", "N Waveform", 11, -0.5, 10.5);
  _hlen = tfs.make<TH1D>("len", "len", 51, -0.5, 50.5);
  _hadc = tfs.make<TH1D>("ADC", "ADC", 100, -10000.0, 2000.0);
  _hamax = tfs.make<TH1D>("ADCmax", "ADCmax", 100, 0.0, 10000.0);

  return 0;
}

int mu2e::ValSTMWaveformDigi::fill(
    const mu2e::STMWaveformDigiCollection& coll, art::Event const& event) {
  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hNwf->Fill(coll.size());
  for (auto wf : coll) {
    auto const& adcs = wf.adcs();
    _hlen->Fill(adcs.size());
    int max = 0;
    for (auto const a : adcs) {
      max = std::max(max, a - adcs[0]);
      _hadc->Fill(a);
    }
    _hamax->Fill(max);  // ADC value as peak minus pedestal
  }
  return 0;
}
