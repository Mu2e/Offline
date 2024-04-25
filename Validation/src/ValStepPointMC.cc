
#include "Offline/Validation/inc/ValStepPointMC.hh"

int mu2e::ValStepPointMC::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>("NStep", "N Steps", 100, -0.05, 1000.0);
  _id.declare(tfs);
  _hp = tfs.make<TH1D>("p", "P", 100, 0.0, 200.0);
  _het1 = tfs.make<TH1D>("eTot1", "E total", 100, 0.0, 2.00);
  _het2 = tfs.make<TH1D>("eTot2", "E total", 100, 0.0, 0.010);
  _ht = tfs.make<TH1D>("t", "time", 100, 0.0, 2000.0);
  _hx = tfs.make<TH1D>("X", "X", 100, -6000.0, 6000.0);
  _hy = tfs.make<TH1D>("Y", "Y", 100, -2000.0, 2000.0);
  _hz = tfs.make<TH1D>("Z", "Z", 100, -20000.0, 20000.0);
  _hl1 = tfs.make<TH1D>("Length1", "Length", 100, 0.0, 100.0);
  _hl2 = tfs.make<TH1D>("Length2", "Length", 100, 0.0, 10.0);
  _hl3 = tfs.make<TH1D>("Length3", "log10(Length)", 100, -19.0, 6.0);
  _hxDS = tfs.make<TH1D>("XDS", "X DS", 100, -3904 - 1000, -3904 + 1000);
  _hyDS = tfs.make<TH1D>("YDS", "Y DS", 100, -1000.0, 1000.0);
  _hxTrk = tfs.make<TH1D>("XTrk", "X Trk", 100, -1000.0, 1000.0);
  _hzTrk = tfs.make<TH1D>("ZTrk", "Z Trk", 100, -1600.0, 1600.0);
  _hzCal = tfs.make<TH1D>("ZCal", "Z Cal", 100, 11750.0, 12900.0);
  _hxCRV = tfs.make<TH1D>("XCRV", "X CRV", 100, -7500.0, 0.0);
  _hyCRV = tfs.make<TH1D>("YCRV", "Y CRV", 100, -2000.0, 3000.0);
  _hzCRV = tfs.make<TH1D>("ZCRV", "Z CRV", 100, -4000.0, 20000.0);

  return 0;
}

int mu2e::ValStepPointMC::fill(const mu2e::StepPointMCCollection& coll,
                               art::Event const& event) {
  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(1.0);

  _hN->Fill(coll.size());

  for (auto sp : coll) {
    // unfortunately, if the event was truncated, then StepPointMC's
    // can point to missing SimParticles, this code catches that case
    // (none of the art::Ptr checks will tell us this, those only abort)
    if (sp.simParticle().isAvailable()) {
      _id.fill(sp.simParticle()->pdgId());
    }

    _hp->Fill(sp.momentum().mag());
    _hx->Fill(sp.position().x());
    _hy->Fill(sp.position().y());
    _hz->Fill(sp.position().z());
    _hl1->Fill(sp.stepLength());
    _hl2->Fill(sp.stepLength());
    _hl3->Fill(log10(sp.stepLength()));
    _hxDS->Fill(sp.position().x());
    _hyDS->Fill(sp.position().y());
    _hxTrk->Fill(sp.position().x());
    _hzTrk->Fill(sp.position().z());
    _hzCal->Fill(sp.position().z());
    _hxCRV->Fill(sp.position().x());
    _hyCRV->Fill(sp.position().y());
    _hzCRV->Fill(sp.position().z());
    _het1->Fill(sp.totalEDep());
    _het2->Fill(sp.totalEDep());
    _ht->Fill(sp.time());
  }
  return 0;
}
