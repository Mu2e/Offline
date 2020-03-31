
#include "Validation/inc/ValComboHit.hh"

int mu2e::ValComboHit::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NHit", "N Combo Hits", 101, -0.5, 100.0);
  _hN2 = tfs.make<TH1D>( "NHit2", "N Combo Hits", 100, -0.5, 9999.5);
  _hNstr = tfs.make<TH1D>( "Nstr", "N Parent Straws", 11, -0.5, 10.5);
  _hWD = tfs.make<TH1D>( "WD", "WireDistance;Distance from mid (mm)", 100, -800.0, 800.0);
  _hDE = tfs.make<TH1D>( "DE", "TOT Drift Time Estimate;Drift Time (ns)", 100, 0.0, 80.0);
  _ht = tfs.make<TH1D>( "t", "time;Hit Time (ns)", 100, 0.0, 2000.0);
  _hE = tfs.make<TH1D>( "E", "Energy;Hit Energy (KeV)",50, 0.0, 5.0);
  _hqual = tfs.make<TH1D>( "qual", "Qual",50, 0.001, 3.0);
  _hwres = tfs.make<TH1D>( "wres", "Wire resolution;resolution (mm)",50, 0.0, 100.0);
  _htres = tfs.make<TH1D>( "tres", "Transverse resolution;resolution (mm)",50, 0.0, 10.0);
  _hPanel = tfs.make<TH1D>( "Panel", "Tracker Unique Panel",216, -0.5, 215.5);
  _hStraw = tfs.make<TH1D>( "Straw", "Tracker Straw",96, -0.5, 95.5);

  return 0;
}

int mu2e::ValComboHit::fill(const mu2e::ComboHitCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(1.0);

  _hN->Fill(coll.size()); 
  _hN2->Fill(coll.size());
  for(auto co : coll) {
    _hNstr->Fill(co.nStrawHits());
    _hWD->Fill(co.wireDist());
    _hDE->Fill(co.driftTime());
    _ht->Fill(co.time());
    _hE->Fill(1000*co.energyDep());
    _hqual->Fill(co.qual());
    _hwres->Fill(co.wireRes());
    _htres->Fill(co.transRes());
    _hPanel->Fill(co.strawId().uniquePanel());
    _hStraw->Fill(co.strawId().straw());
  }
  return 0;
}
