
#include "Validation/inc/ValStereoHit.hh"

int mu2e::ValStereoHit::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NHit", "N Stereo Hits", 101, -0.5, 100.0);
  _hN2 = tfs.make<TH1D>( "NHit2", "N Stereo Hits", 100, -0.5, 9999.5);
  _pos.declare(tfs,"","");
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _hE = tfs.make<TH1D>( "E", "Energy",50, 0.0, 0.02);
  _hchi2 = tfs.make<TH1D>( "chi2", "Chi2",50, 0.0, 50.0);
  _hmva = tfs.make<TH1D>( "mva", "MVA",50, 0.0, 1.00);

  return 0;
}

int mu2e::ValStereoHit::fill(const mu2e::StereoHitCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size()); 
  _hN2->Fill(coll.size()); 
  for(auto sh : coll) {
    _pos.fill(sh.pos());
    _ht->Fill(sh.time());
    _hE->Fill(sh.energy());
    _hchi2->Fill(sh.chisq());
    _hmva->Fill(sh.mvaout());
  }
  return 0;
}
