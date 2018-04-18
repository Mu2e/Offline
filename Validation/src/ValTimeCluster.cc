
#include "Validation/inc/ValTimeCluster.hh"

int mu2e::ValTimeCluster::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "N", "N Time Clusters", 101, -0.5, 100.0);
  _hNhit = tfs.make<TH1D>( "NHit", "N Hits", 101, -0.5, 100.0);
  _pos.declare(tfs,"","");
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);

  return 0;
}

int mu2e::ValTimeCluster::fill(const mu2e::TimeClusterCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size()); 
  for(auto tp : coll) {
    _hNhit->Fill(tp.nhits());
    CLHEP::Hep3Vector pos = Geom::Hep3Vec(tp.position());
    _pos.fill(pos);
    _ht->Fill(tp.t0().t0());
  }
  return 0;
}
