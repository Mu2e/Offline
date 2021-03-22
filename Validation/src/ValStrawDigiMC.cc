
#include "Validation/inc/ValStrawDigiMC.hh"

int mu2e::ValStrawDigiMC::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NHit", "N Straw Hits", 101, -0.5, 100.0);
  _hN2 = tfs.make<TH1D>( "NHit2", "N Straw Hits", 100, -0.5, 9999.5);
  _htime0 = tfs.make<TH1D>( "time0", "time0", 100, 0.0, 2000.0);
  _htime1 = tfs.make<TH1D>( "time1", "time1", 100, 0.0, 2000.0);
  _hener = tfs.make<TH1D>( "ener", "energy",50, -1.0e-6, 0.01);
  _henerT = tfs.make<TH1D>( "enerT", "trigger energy",50, -1.0e-6, 0.01);
  _hcross = tfs.make<TH1D>( "cross", "cross talk flag",2, -0.5, 1.5);
  _hgStep = tfs.make<TH1D>( "gStep", "N good step pointers",3, -0.5, 2.5);
  _hSI = tfs.make<TH1D>( "Straw", "Straw Index",100, -0.5, StrawId::_maxval);

  return 0;
}

int mu2e::ValStrawDigiMC::fill(const mu2e::StrawDigiMCCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(3.0);

  _hN->Fill(coll.size()); 
  _hN2->Fill(coll.size()); 


  for(auto sd : coll) {

    // check the Ptr's.  If the products are not there, the accessors can crash
    bool ptrOK = true;
    auto const& a0 = sd.strawGasStep(StrawEnd::cal);
    auto const& a1 = sd.strawGasStep(StrawEnd::hv);
    
    if(!(a0.isAvailable() && a1.isAvailable()) ) ptrOK = false;

    double ns = 0.0;
    if(sd.strawGasStep(StrawEnd::cal).isAvailable()) ns++;
    if(sd.strawGasStep(StrawEnd::hv).isAvailable()) ns++;
    _hgStep->Fill(ns);
    
    if(ptrOK) {
        _htime0->Fill(sd.wireEndTime(StrawEnd::cal));
      _htime1->Fill(sd.wireEndTime(StrawEnd::hv));
      _hener->Fill(sd.energySum());
      _henerT->Fill(sd.triggerEnergySum(StrawEnd::cal));
      _henerT->Fill(sd.triggerEnergySum(StrawEnd::hv));
      if(sd.isCrossTalk(StrawEnd::cal)) _hcross->Fill(0.0);
      if(sd.isCrossTalk(StrawEnd::hv) ) _hcross->Fill(1.0);
      _hSI->Fill(sd.strawId().asUint16()); 
    }
  }
  return 0;
}
