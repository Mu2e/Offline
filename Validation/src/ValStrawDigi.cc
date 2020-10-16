
#include "Validation/inc/ValStrawDigi.hh"
#include "DataProducts/inc/StrawEnd.hh"

int mu2e::ValStrawDigi::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NHit", "N Straw Hits", 101, -0.5, 100.0);
  _hN2 = tfs.make<TH1D>( "NHit2", "N Straw Hits", 100, -0.5, 9999.5);
  _htdc = tfs.make<TH1D>( "TDC", "TDC", 100, 0.0, 120000.0);
  _hadc = tfs.make<TH1D>( "ADC", "ADC",50, -0.5, 5000.0);
  _hSI = tfs.make<TH1D>( "Straw", "Straw Index",100, -0.5, 40960.0);

  return 0;
}

int mu2e::ValStrawDigi::fill(const mu2e::StrawDigiCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(1.0);

  _hN->Fill(coll.size()); 
  _hN2->Fill(coll.size()); 
  for(auto sd : coll) {
    for (size_t ie=0; ie<StrawEnd::nends; ie++) {  // for the two straw ends
      _htdc->Fill(sd.TDC()[ie]);
    }
    //for(auto const& a : sd.adcWaveform()) _hadc->Fill(a); // ADC values (ints)
    //_hSI->Fill(sd.strawId().asUint16()); // <40960
  }
  return 0;
}
