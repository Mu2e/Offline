
#include "Validation/inc/ValTimeCluster.hh"

int mu2e::ValTimeCluster::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "N", "N Time Clusters", 101, -0.5, 100.0);
  _hNhit = tfs.make<TH1D>( "NHit", "N Hits", 101, -0.5, 100.0);
  _hx = tfs.make<TH1D>( "X", "tracker X", 100, -800.0, 800.0);
  _hy = tfs.make<TH1D>( "Y", "tracker Y", 100, -800.0, 800.0);
  _hz = tfs.make<TH1D>( "Z", "tracker Z", 100, -1500.0, 3000.0);
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
    _hx->Fill(pos.x());
    _hy->Fill(pos.y());
    _hz->Fill(pos.z());
    _ht->Fill(tp.t0().t0());
  }
  return 0;
}
