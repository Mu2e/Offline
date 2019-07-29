
#include "Validation/inc/ValBkgQual.hh"


int mu2e::ValBkgQual::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NClus", "N Clusters", 100, 0.001, 500.0);
  _hmva = tfs.make<TH1D>( "mva", "mva", 50, -0.02, 1.02);
  _hstat = tfs.make<TH1D>( "bits", "BkgClusterFlag", 11, -0.5, 10.5);

  return 0;
}

int mu2e::ValBkgQual::fill(const mu2e::BkgQualCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);
  _hN->Fill(coll.size()); 
  for(auto bc : coll) {
    _hmva->Fill(bc.MVAOutput());
    _hstat->Fill((int)bc.status());
  }
  return 0;
}
