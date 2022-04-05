
#include "Offline/Validation/inc/ValProtonBunchTimeMC.hh"


int mu2e::ValProtonBunchTimeMC::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>( "Ver", "Version Number", 101, -0.5, 100.0);
  _htime = tfs.make<TH1D>( "time", "time", 50, -250.0, -150.0);

  return 0;
}

int mu2e::ValProtonBunchTimeMC::fill(const mu2e::ProtonBunchTimeMC & obj,
                                art::Event const& event) {

  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _htime->Fill(obj.pbtime_);

  return 0;
}
