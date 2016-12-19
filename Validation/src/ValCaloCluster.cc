
#include "Validation/inc/ValCaloCluster.hh"


int mu2e::ValCaloCluster::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NClus", "N Clusters", 31, -0.05, 30.0);
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _hE = tfs.make<TH1D>( "E", "Energy",50, 0.0, 200.0);
  _hR = tfs.make<TH1D>( "R", "Radius",50, 340.0, 650.0);
  _hA = tfs.make<TH1D>( "A", "angle", 50, 0.0, 3.1416);

  return 0;
}

int mu2e::ValCaloCluster::fill(const mu2e::CaloClusterCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

   _hN->Fill(coll.size()); 
  for(auto sp : coll) {
    _ht->Fill(sp.time());
    _hE->Fill(sp.energyDep());
    _hR->Fill(sp.cog3Vector().perp());
    _hA->Fill(sp.angle());
  }
  return 0;
}
