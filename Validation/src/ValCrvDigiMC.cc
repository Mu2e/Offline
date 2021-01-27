
#include "Validation/inc/ValCrvDigiMC.hh"


int mu2e::ValCrvDigiMC::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NDigis", "N Digis", 101, -0.5, 100.5);
  _hN2= tfs.make<TH1D>( "NDigis2", "N Digis", 100, -0.5, 4999.5);
  _hI = tfs.make<TH1D>( "BarId", "Bar ID",200, -0.5, 5503.5);
  _hIS= tfs.make<TH1D>( "SiPM", "SiPM",4, -0.5, 3.5);
  _hNS= tfs.make<TH1D>( "NStep", "N Step Points",31, -0.5, 30.5);
  _ht = tfs.make<TH1D>( "t", "time", 100, 0.0, 2000.0);
  _hV = tfs.make<TH1D>( "V", "Voltages",100, 0.0, 6.0);

  return 0;
}

int mu2e::ValCrvDigiMC::fill(const mu2e::CrvDigiMCCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size());
  _hN2->Fill(coll.size());
  for(auto dg : coll) {
    _hI->Fill(dg.GetScintillatorBarIndex().asInt());
    _hIS->Fill(dg.GetSiPMNumber());
    _hNS->Fill(dg.GetCrvSteps().size());
    _ht->Fill(dg.GetStartTime());
    _ht->Fill(dg.GetStartTime());
    for (auto v: dg.GetVoltages()) _hV->Fill(v);

  }
  return 0;
}
