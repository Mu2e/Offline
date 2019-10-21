
#include "Validation/inc/ValTriggerResults.hh"


int mu2e::ValTriggerResults::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hNpath = tfs.make<TH1D>( "Npath", "N path", 50, -0.5, 50.5);
  _hState = tfs.make<TH1D>( "State", "PathNo*4+State(i)", 201, -0.5, 200.5);
  _hIndex = tfs.make<TH1D>( "Index", "Index for all paths", 101, -0.5, 100.5);

  return 0;
}

int mu2e::ValTriggerResults::fill(const art::TriggerResults & obj,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hNpath->Fill(obj.size());
  for(size_t i=0; i<obj.size(); i++) {
    _hState->Fill(i*4+obj.state(i));
    _hIndex->Fill(obj.index(i));
  }

  return 0;
}
