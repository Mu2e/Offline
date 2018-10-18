
#include "Validation/inc/ValCrvRecoPulse.hh"


int mu2e::ValCrvRecoPulse::declare(art::TFileDirectory tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _hN = tfs.make<TH1D>( "NPulses", "N Pulses", 101, -0.5, 100.5);
  _hN2= tfs.make<TH1D>( "NPulse2", "N Pulses", 100, -0.5, 2000.0);
  _hI = tfs.make<TH1D>( "BarId", "Bar ID",200, -0.5, 5503.5);
  _hIS= tfs.make<TH1D>( "SiPM", "SiPM",4, -0.5, 3.5);
  _hPE = tfs.make<TH1D>( "PE", "Fit Photoelectrons",100, 0.0, 400.0);
  _hPH = tfs.make<TH1D>( "PEHeight", "PE from Pulse Height", 100, 0.0, 400.0);
  _ht  = tfs.make<TH1D>( "PulseTime", "Pulse Peak Time", 100, 0.0, 2000.0);
  _hChi2 = tfs.make<TH1D>( "chi2", "Pulse fit chi2",100, 0.0, 1.0e5);
  _hLE = tfs.make<TH1D>( "LeadingTime", "Leading Edge Time",100, 0.0, 2000.0);

  return 0;
}

int mu2e::ValCrvRecoPulse::fill(const mu2e::CrvRecoPulseCollection & coll,
				art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the 
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hN->Fill(coll.size());
  _hN2->Fill(coll.size());
  for(auto rp : coll) {
    _hI->Fill(rp.GetScintillatorBarIndex().asInt());
    _hIS->Fill(rp.GetSiPMNumber());
    _hPE->Fill(rp.GetPEs());
    _hPH->Fill(rp.GetPEsPulseHeight());
    _ht->Fill(rp.GetPulseTime());
    _hChi2->Fill(rp.GetPulseFitChi2());
    _hLE->Fill(rp.GetLEtime());

  }

  return 0;
}
