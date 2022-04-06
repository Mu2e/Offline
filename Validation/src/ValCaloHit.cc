
#include "Offline/Validation/inc/ValCaloHit.hh"

int mu2e::ValCaloHit::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>("NDigis", "N Hits", 101, -0.5, 100.5);
  _hN2 = tfs.make<TH1D>("NDigis2", "N Hits", 100, -0.5, 4999.5);
  _hI = tfs.make<TH1D>("crystalID", "crystalID", 100, 0.0, 1400.0);
  _ht = tfs.make<TH1D>("t", "time", 100, 0.0, 2000.0);
  _hE = tfs.make<TH1D>("E", "Energy", 100, 0.0, 150.0);

  return 0;
}

int mu2e::ValCaloHit::fill(const mu2e::CaloHitCollection& coll,
                           art::Event const& event) {
  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size());
  _hN2->Fill(coll.size());
  for (auto ch : coll) {
    _hI->Fill(ch.crystalID());
    _ht->Fill(ch.time());
    _hE->Fill(ch.energyDep());
  }
  return 0;
}
