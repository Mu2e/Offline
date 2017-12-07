
#include "Validation/inc/ValStepPointMC.hh"


int mu2e::ValStepPointMC::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NStep", "N Steps", 100, -0.05, 1000.0);
  _id.declare(tfs);
  _hp = tfs.make<TH1D>( "p", "P", 100, 0.0, 200.0);
  _het1 = tfs.make<TH1D>( "eTot1", "E total", 100, 0.0, 2.00);
  _het2 = tfs.make<TH1D>( "eTot2", "E total", 100, 0.0, 0.010);
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _pos.declare(tfs);
  return 0;
}

int mu2e::ValStepPointMC::fill(const mu2e::StepPointMCCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size()); 
  for(auto sp : coll) {
    _id.fill(sp.simParticle()->pdgId()); 
    _hp->Fill(sp.momentum().mag()); 
    _pos.fill(sp.position());
    _het1->Fill(sp.totalEDep()); 
    _het2->Fill(sp.totalEDep()); 
    _ht->Fill(sp.time()); 
  }
  return 0;
}
