#include "Offline/Validation/inc/ValStrawGasStep.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

int mu2e::ValStrawGasStep::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NHit", "N Straw Hits", 101, -0.5, 100.0);
  _hN2 = tfs.make<TH1D>( "NHit2", "N Straw Hits", 100, -0.5, 9999.5);
  _ht = tfs.make<TH1D>( "t", "time", 100, -100.0, 2000.0);
  _ht2 = tfs.make<TH1D>( "t2", "time", 100, -100.0, 100.0e3);
  _hE = tfs.make<TH1D>( "E", "Energy",50, 0.0, 0.01);
  _hlE = tfs.make<TH1D>( "lE", "log10(Energy)",100, -6.0, 1.0);
  _hlen = tfs.make<TH1D>( "Length", "steplength",100, 0.0, 10.0);
  _pmom = tfs.make<TH1D>( "log10(pmom)", "Particle momentum",100, -5,5);
  _hz = tfs.make<TH1D>( "Z", "Z",100, -1600.0, 1600.0);
  _hSI = tfs.make<TH1D>( "Straw", "Unique Straw",100, -0.5, StrawId::_nustraws-0.5);
  _hpla = tfs.make<TH1D>( "Plane", "Plane",36, -0.5, StrawId::_nplanes-0.5);
  _hstr = tfs.make<TH1D>( "RStraw", "Reduced straw",96, -0.5, StrawId::_nstraws-0.5);

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
    _ht2->Fill(gs.time());
    _hE->Fill(gs.ionizingEdep());
    _hlE->Fill(  ( gs.ionizingEdep()>0.0 ? log10(gs.ionizingEdep()) : -10 ) );
    _hlen->Fill(gs.stepLength());
    _pmom->Fill(log10(gs.momentum().R()));
    _hz->Fill(gs.startPosition().z());
    _hSI->Fill(gs.strawId().uniqueStraw());
    _hpla->Fill(gs.strawId().plane());
    _hstr->Fill(gs.strawId().straw());
  }
  return 0;
}
