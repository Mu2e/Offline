
#include "Offline/Validation/inc/ValStrawDigiADCWaveform.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"
#include <algorithm>

int mu2e::ValStrawDigiADCWaveform::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>("NHit", "N Waveform", 101, -0.5, 100.0);
  _hN2 = tfs.make<TH1D>("NHit2", "N Waveform", 100, -0.5, 9999.5);
  _hlen = tfs.make<TH1D>("len", "len", 51, -0.5, 50.5);
  _hadc = tfs.make<TH1D>("ADC", "ADC", 100, 200.0, 1200.0);
  _hpmp = tfs.make<TH1D>("PMP", "PMP", 100, -0.5, 1000.0);

  return 0;
}

int mu2e::ValStrawDigiADCWaveform::fill(
    const mu2e::StrawDigiADCWaveformCollection& coll, art::Event const& event) {
  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size());
  _hN2->Fill(coll.size());
  for (auto wf : coll) {
    _hlen->Fill(wf.samples().size());
    int max = 0;
    for (auto const a : wf.samples()) {  // for the two straw ends
      max = std::max(max, int(a) - int(wf.samples()[0]));
      _hadc->Fill(a);
    }
    _hpmp->Fill(max);  // ADC value as peak minus pedestal
  }
  return 0;
}
