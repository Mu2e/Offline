
#include "Validation/inc/ValTimeCluster.hh"

int mu2e::ValTimeCluster::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.5);
  _hN = tfs.make<TH1D>( "N", "N Time Clusters", 101, -0.5, 100.5);
  _hNhit = tfs.make<TH1D>( "NHit", "N Hits", 101, -0.5, 100.5);
  _hx = tfs.make<TH1D>( "X", "tracker X", 100, -800.0, 800.0);
  _hy = tfs.make<TH1D>( "Y", "tracker Y", 100, -800.0, 800.0);
  _hz = tfs.make<TH1D>( "Z", "tracker Z", 100, -1500.0, 3000.0);
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _hterr = tfs.make<TH1D>( "terr", "time error", 100, 0.0, 20.0);
  _nc = tfs.make<TH1D>( "nc", "N Calo Clusters", 2, -0.5, 1.5);

  return 0;
}

int mu2e::ValTimeCluster::fill(const mu2e::TimeClusterCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(1.0);

  _hN->Fill(coll.size()); 
  for(auto tp : coll) {
    _hNhit->Fill(tp.nStrawHits());
    CLHEP::Hep3Vector pos = Geom::Hep3Vec(tp.position());
    _hx->Fill(pos.x());
    _hy->Fill(pos.y());
    _hz->Fill(pos.z());
    _ht->Fill(tp.t0().t0());
    _hterr->Fill(tp.t0().t0Err());
    if(tp.hasCaloCluster())
      _nc->Fill(1.0);
    else
      _nc->Fill(0.0);
  }
  return 0;
}
