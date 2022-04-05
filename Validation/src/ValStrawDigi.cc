
#include "Offline/Validation/inc/ValStrawDigi.hh"
#include "Offline/DataProducts/inc/StrawEnd.hh"

int mu2e::ValStrawDigi::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NHit", "N Straw Hits", 101, -0.5, 100.0);
  _hN2 = tfs.make<TH1D>( "NHit2", "N Straw Hits", 100, -0.5, 9999.5);
  _htdc = tfs.make<TH1D>( "TDC", "TDC", 100, 0.0, 80.0e3);
  _htdc2= tfs.make<TH1D>( "TDC2", "TDC2", 100, 0.0, 5.0e6);
  _htot = tfs.make<TH1D>( "TOT", "TOT", 21, -0.5, 20.5);
  _hpmp = tfs.make<TH1D>( "PMP", "PMP",50, -0.5, 1000.0);
  _hSI = tfs.make<TH1D>( "Straw", "Straw Index",100, -0.5, StrawId::_maxval);

  return 0;
}

int mu2e::ValStrawDigi::fill(const mu2e::StrawDigiCollection & coll,
                                art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(2.0);

  _hN->Fill(coll.size());
  _hN2->Fill(coll.size());
  for(auto sd : coll) {
    for (size_t ie=0; ie<StrawEnd::nends; ie++) {  // for the two straw ends
      _htdc->Fill(sd.TDC()[ie]);
      _htdc2->Fill(sd.TDC()[ie]);
      _htot->Fill(sd.TOT()[ie]);
    }
    _hpmp->Fill(sd.PMP()); // ADC value as peak minus pedestal
    _hSI->Fill(sd.strawId().asUint16());
  }
  return 0;
}
