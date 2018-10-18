
#include "Validation/inc/ValComboHit.hh"

int mu2e::ValComboHit::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NHit", "N Stereo Hits", 101, -0.5, 100.0);
  _hN2 = tfs.make<TH1D>( "NHit2", "N Stereo Hits", 100, -0.5, 9999.5);
  _hsid = tfs.make<TH1D>( "sid", "StrawId", 100, -0.5, 40960.0);
  _hNcmb = tfs.make<TH1D>( "Ncmb", "N Parent Combos", 11, -0.5, 10.5);
  _hNstr = tfs.make<TH1D>( "Nstr", "N Parent Straws", 11, -0.5, 10.5);
  _hx = tfs.make<TH1D>( "X", "tracker X", 100, -800.0, 800.0);
  _hy = tfs.make<TH1D>( "Y", "tracker Y", 100, -800.0, 800.0);
  _hz = tfs.make<TH1D>( "Z", "tracker Z", 100, -1600.0, 1600.0);
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _hE = tfs.make<TH1D>( "E", "Energy",50, 0.0, 0.02);
  _hqual = tfs.make<TH1D>( "qual", "Qual",50, 0.001, 3.0);
  _hwres = tfs.make<TH1D>( "wres", "Wire resolution",50, 0.0, 100.0);
  _htres = tfs.make<TH1D>( "tres", "Transvrse resolution",50, 0.0, 10.0);

  return 0;
}

int mu2e::ValComboHit::fill(const mu2e::ComboHitCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size()); 
  _hN2->Fill(coll.size());
  for(auto co : coll) {
    _hsid->Fill(co.strawId().asUint16());
    _hNcmb->Fill(co.nCombo());
    _hNstr->Fill(co.nStrawHits());
    _hx->Fill(co.posCLHEP().x());
    _hy->Fill(co.posCLHEP().y());
    _hz->Fill(co.posCLHEP().z());
    _ht->Fill(co.time());
    _hE->Fill(co.energyDep());
    _hqual->Fill(co.qual());
    _hwres->Fill(co.wireRes());
    _htres->Fill(co.transRes());
  }
  return 0;
}
