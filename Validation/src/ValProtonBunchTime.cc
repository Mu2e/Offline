
#include "Offline/Validation/inc/ValProtonBunchTime.hh"

int mu2e::ValProtonBunchTime::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _htime = tfs.make<TH1D>("time", "time", 50, -250.0, -150.0);
  _hterr = tfs.make<TH1D>("terr", "time error", 50, 0.0, 20.0);

  return 0;
}

int mu2e::ValProtonBunchTime::fill(const mu2e::ProtonBunchTime& obj,
                                   art::Event const& event) {
  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(1.0);

  _htime->Fill(obj.pbtime_);
  _hterr->Fill(obj.pbterr_);

  return 0;
}
