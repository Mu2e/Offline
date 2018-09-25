
#include "Validation/inc/ValCrvDigi.hh"


int mu2e::ValCrvDigi::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NDigis", "N Digis", 101, -0.5, 100.5);
  _hN2= tfs.make<TH1D>( "NDigis2", "N Digis", 100, -0.5, 4999.5);
  _hI = tfs.make<TH1D>( "BarId", "Bar ID",200, -0.5, 5503.5);
  _hIS= tfs.make<TH1D>( "SiPM", "SiPM",4, -0.5, 3.5);
  _ht = tfs.make<TH1D>( "t", "TDC", 100, 0.0, 250.0);
  _hA = tfs.make<TH1D>( "ADC", "ADC in waveform",100, 0.0, 3000.0);

  return 0;
}

int mu2e::ValCrvDigi::fill(const mu2e::CrvDigiCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size());
  _hN2->Fill(coll.size());
  for(auto dg : coll) {
    _hI->Fill(dg.GetScintillatorBarIndex().asInt());
    _hIS->Fill(dg.GetSiPMNumber());
    _ht->Fill(dg.GetStartTDC());
    for (auto a: dg.GetADCs()) _hA->Fill(a);
  }
  return 0;
}
