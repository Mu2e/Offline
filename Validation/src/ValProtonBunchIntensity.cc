
#include "Offline/Validation/inc/ValProtonBunchIntensity.hh"

int mu2e::ValProtonBunchIntensity::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _hint = tfs.make<TH1D>("int", "PBI/1e6", 100, 0.0, 200.0);

  return 0;
}

int mu2e::ValProtonBunchIntensity::fill(const mu2e::ProtonBunchIntensity& obj,
                                        art::Event const& event) {
  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);

  _hint->Fill(obj.intensity() / 1.0e6);

  return 0;
}
