
#include "Validation/inc/ValCaloCrystalHit.hh"


int mu2e::ValCaloCrystalHit::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NHits", "N Hits", 101, -0.5, 100.5);
  _hI = tfs.make<TH1D>( "ID", "ID",150, 0.0, 1500.0);
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _hE = tfs.make<TH1D>( "E", "Energy",100, 0.0, 200.0);

  return 0;
}

int mu2e::ValCaloCrystalHit::fill(const mu2e::CaloCrystalHitCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

   _hN->Fill(coll.size()); 
  for(auto sp : coll) {
    _hI->Fill(sp.id());
    _ht->Fill(sp.time());
    _hE->Fill(sp.energyDep());
  }
  return 0;
}
