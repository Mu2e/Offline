
#include "Validation/inc/ValCaloShowerStep.hh"


int mu2e::ValCaloShowerStep::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NStep", "N Steps", 101, -9.5, 1000.5);
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _hE = tfs.make<TH1D>( "E", "Energy",50, 0.0, 200.0);
  _hposx = tfs.make<TH1D>( "posx", "Position X",50, -20.0, 20.0);
  _hposy = tfs.make<TH1D>( "posy", "Position Y",50, -20.0, 20.0);
  _hposz = tfs.make<TH1D>( "posz", "Position Z",50, 0.0, 210.0);


  return 0;
}

int mu2e::ValCaloShowerStep::fill(const mu2e::CaloShowerStepCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size()); 
  for(auto ss : coll) {
    _ht->Fill(ss.timeStepMC());
    _hE->Fill(ss.energyMC());
    _hposx->Fill(ss.position().x());
    _hposy->Fill(ss.position().y());
    _hposz->Fill(ss.position().z());
  }
  return 0;
}
