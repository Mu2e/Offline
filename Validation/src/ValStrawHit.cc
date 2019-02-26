
#include "Validation/inc/ValStrawHit.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"

int mu2e::ValStrawHit::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NHit", "N Straw Hits", 101, -0.5, 100.0);
  _hN2 = tfs.make<TH1D>( "NHit2", "N Straw Hits", 100, -0.5, 9999.5);
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _hdt = tfs.make<TH1D>( "dt", "d-time", 100, -15.0, 15.0);
  _hE = tfs.make<TH1D>( "E", "Energy",50, 0.0, 0.01);
  _hDI = tfs.make<TH1D>( "PlanePanel", "Plane*6+Panel Index",240, -0.5, 239.5);
  _hSI = tfs.make<TH1D>( "Straw", "Straw Index",96, -0.5, 95.5);

  return 0;
}

int mu2e::ValStrawHit::fill(const mu2e::StrawHitCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  const mu2e::Tracker& tracker = *GeomHandle<mu2e::Tracker>();

  _hN->Fill(coll.size()); 
  _hN2->Fill(coll.size()); 
  for(auto sp : coll) {
    Straw const& straw = tracker.getStraw( sp.strawId() );
    StrawId const& id = straw.id();
    _ht->Fill(sp.time());
    _hdt->Fill(sp.dt());
    _hE->Fill(sp.energyDep());
    _hDI->Fill(6*id.getPlane()+id.getPanel());
    _hSI->Fill(id.getStraw());
  }
  return 0;
}
