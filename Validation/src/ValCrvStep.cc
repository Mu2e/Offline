
#include "Offline/Validation/inc/ValCrvStep.hh"

int mu2e::ValCrvStep::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>("NStep", "N Steps", 101, -9.5, 1000.5);
  _hb = tfs.make<TH1D>("bar", "bar number", 100, -0.5, 5503.5);
  _ht = tfs.make<TH1D>("t", "time", 100, 0.0, 2000.0);
  _ht2 = tfs.make<TH1D>("t2", "time", 100, -100.0, 100.0e3);
  _hlt =
      tfs.make<TH1D>("lt", "log_{10} time;log_{10} time (ns)", 100, 0.0, 10.0);
  _hE = tfs.make<TH1D>("E", "Energy", 100, 0.0, 0.2);
  _hlE = tfs.make<TH1D>("lE", "log10(Energy)", 100, -9.0, 2.0);
  _hposx = tfs.make<TH1D>("posx", "Start position X", 100, -7500.0, 1000.0);
  _hposy = tfs.make<TH1D>("posy", "Start position Y", 100, -2000.0, 3000.0);
  _hposz = tfs.make<TH1D>("posz", "Start position Z", 100, -4000.0, 20000.0);
  _hp = tfs.make<TH1D>("SMom", "Start Momentum", 100, 0.0, 100.0);
  _hp2 = tfs.make<TH1D>("SMom2", "Start Momentum 2", 100, 0.0, 100000.0);
  _hpL = tfs.make<TH1D>("pathLen", "pathLength", 100, 0.0, 100.0);
  return 0;
}

int mu2e::ValCrvStep::fill(const mu2e::CrvStepCollection& coll,
                           art::Event const& event) {
  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(2.0);

  _hN->Fill(coll.size());
  for (auto cs : coll) {
    _hb->Fill(cs.barIndex().asInt());
    _ht->Fill(cs.startTime());
    _ht2->Fill(cs.startTime());
    _hlt->Fill(log10(std::max(cs.startTime(), 1e-5)));
    _hE->Fill(cs.visibleEDep());
    _hlE->Fill((cs.visibleEDep() > 0.0 ? log10(cs.visibleEDep()) : -10.0));
    _hposx->Fill(cs.startPosition().x());
    _hposy->Fill(cs.startPosition().y());
    _hposz->Fill(cs.startPosition().z());
    _hp->Fill(cs.startMomentum().mag());
    _hp2->Fill(cs.startMomentum().mag());
    _hpL->Fill(cs.pathLength());
  }
  return 0;
}
