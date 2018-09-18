
#include "Validation/inc/ValStrawHitFlag.hh"
#include <cmath>

int mu2e::ValStrawHitFlag::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NFlag", "N Straw Hit flags", 101, -0.5, 100.5);
  _hN2 = tfs.make<TH1D>( "NFlag2", "N Straw Hit flags", 100, -0.5, 9999.5);
  _hBits = tfs.make<TH1D>( "Bits", "Bits", 32, -0.5, 31.5);

  return 0;
}

int mu2e::ValStrawHitFlag::fill(const mu2e::StrawHitFlagCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size()); 
  _hN2->Fill(coll.size()); 
  for(auto f : coll) {
    for(auto sn: f.bitNames()) { 
      if(f.hasAnyProperty(StrawHitFlag(sn.first))) _hBits->Fill(std::log2(sn.second)); 
    }
  }
  return 0;
}
