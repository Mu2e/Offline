
#include "Validation/inc/ValBkgCluster.hh"


int mu2e::ValBkgCluster::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NClus", "N Clusters", 100, 0.001, 500.0);
  _hr = tfs.make<TH1D>( "R", "Radius",100, 370.0, 680.0);
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _hd = tfs.make<TH1D>( "d", "hit distance", 100, 0.0, 80.0);
  _hBits = tfs.make<TH1D>( "bits", "BkgClusterFlag", 32, -0.5, 31.5);

  return 0;
}

int mu2e::ValBkgCluster::fill(const mu2e::BkgClusterCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);
  _hN->Fill(coll.size()); 
  for(auto bc : coll) {
    _hr->Fill(bc.pos().R());
    _ht->Fill(bc.time());
    //for(auto const& h : bc.hits()) _hd->Fill(h.distance());    
    for(auto sn: bc.flag().bitNames()) { 
      if(bc.flag().hasAnyProperty(BkgClusterFlag(sn.first))) _hBits->Fill(std::log2(sn.second)); 
    }

  }

  return 0;
}
