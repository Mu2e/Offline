
#include "Offline/Validation/inc/ValCaloDigi.hh"

int mu2e::ValCaloDigi::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>("NDigis", "N Digis", 101, -0.5, 100.5);
  _hN2 = tfs.make<TH1D>("NDigis2", "N Digis", 100, -0.5, 4999.5);
  _hI = tfs.make<TH1D>("SiPMID", "SiPM ID", 100, 0.0, 3000.0);
  _ht = tfs.make<TH1D>("t", "time", 100, 0.0, 2000.0);
  _ht2 = tfs.make<TH1D>("t2", "time", 100, 0.0, 100.0e3);
  _hm = tfs.make<TH1D>("Nwave", "N points in waveform", 81, -0.5, 80.5);
  _hE = tfs.make<TH1D>("EMax", "E max in waveform", 100, 0.0, 5000.0);

  return 0;
}

int mu2e::ValCaloDigi::fill(const mu2e::CaloDigiCollection& coll,
                            art::Event const& event) {
  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size());
  _hN2->Fill(coll.size());
  for (auto dg : coll) {
    _hI->Fill(dg.SiPMID());
    _ht->Fill(dg.t0());
    _ht2->Fill(dg.t0());
    double emax = 0.0;
    _hm->Fill(double(dg.waveform().size()));
    for (auto e : dg.waveform()) {
      if (e > emax) emax = e;
    }
    _hE->Fill(emax);
  }
  return 0;
}
