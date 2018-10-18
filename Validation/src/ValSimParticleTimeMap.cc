
#include "Validation/inc/ValSimParticleTimeMap.hh"

int mu2e::ValSimParticleTimeMap::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "Nsim", "N particle", 101, -0.5, 100.5);
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _ht2 = tfs.make<TH1D>( "t2", "time", 100, -100.0, 100.0);

  return 0;
}

int mu2e::ValSimParticleTimeMap::fill(const mu2e::SimParticleTimeMap & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size()); 
  for(auto sp : coll) {
    _ht->Fill(sp.second);
    _ht2->Fill(sp.second);
  }
  return 0;
}
