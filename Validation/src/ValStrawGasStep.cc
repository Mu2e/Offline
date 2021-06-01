
#include "Validation/inc/ValStrawGasStep.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"

int mu2e::ValStrawGasStep::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NHit", "N Straw Hits", 101, -0.5, 100.0);
  _hN2 = tfs.make<TH1D>( "NHit2", "N Straw Hits", 100, -0.5, 9999.5);
  _ht = tfs.make<TH1D>( "t", "time", 100, -100.0, 2000.0);
  _hE = tfs.make<TH1D>( "E", "Energy",50, 0.0, 0.01);
  _hlen = tfs.make<TH1D>( "Length", "steplength",100, 0.0, 10.0);
  _hz = tfs.make<TH1D>( "Z", "Z",100, -1600.0, 1600.0);
  _hSI = tfs.make<TH1D>( "Straw", "Unique Straw",100, -0.5, StrawId::_nustraws-0.5);

  return 0;
}

int mu2e::ValStrawGasStep::fill(const mu2e::StrawGasStepCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(2.0);

  _hN->Fill(coll.size()); 
  _hN2->Fill(coll.size()); 
  for(auto gs : coll) {
    _ht->Fill(gs.time());
    _hE->Fill(gs.ionizingEdep());
    _hlen->Fill(gs.stepLength());
    _hz->Fill(gs.startPosition().z());
    _hSI->Fill(gs.strawId().uniqueStraw());
  }
  return 0;
}
