
#include "Offline/Validation/inc/ValEventWindowMarker.hh"

int mu2e::ValEventWindowMarker::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _hst = tfs.make<TH1D>("st", "Spill type", 6, -0.5, 5.5);
  _hlen = tfs.make<TH1D>("len", "Event length", 50, 1600.0, 1800.0);

  return 0;
}

int mu2e::ValEventWindowMarker::fill(const mu2e::EventWindowMarker& obj,
                                     art::Event const& event) {
  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hst->Fill(obj.spillType());
  _hlen->Fill(obj.eventLength());

  return 0;
}
