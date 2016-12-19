
#include "Validation/inc/ValCaloHit.hh"


int mu2e::ValCaloHit::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NClus", "N Hits", 151, -0.5, 150.5);
  _hI = tfs.make<TH1D>( "ID", "ID",200, 0.0, 3700.0);
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _hE = tfs.make<TH1D>( "E", "Energy",50, 0.0, 200.0);

  return 0;
}

int mu2e::ValCaloHit::fill(const mu2e::CaloHitCollection & coll,
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
