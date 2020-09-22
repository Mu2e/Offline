#include <memory>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

#include "RecoDataProducts/inc/TrackSummary.hh"
#include "TH1.h"
#include "art_root_io/TFileDirectory.h"

#include "Mu2eUtilities/inc/HistTrackSum.hh"

namespace mu2e {

  HistTrackSum::HistTrackSum(art::TFileDirectory topdir, const std::string& subdir) {
    art::TFileDirectory tfdir = subdir.empty() ? topdir : topdir.mkdir(subdir.c_str());
    TH1::SetDefaultSumw2();
    nactive_ = tfdir.make<TH1D>("nactive", "nactive", 150, -0.5, 149.5);
    fitcon_ = tfdir.make<TH1D>("fitcon", "fitcon", 1000, 0., 1.);
    momerr_ = tfdir.make<TH1D>("momerr", "momerr", 100, 0., 2.);
    t0err_ = tfdir.make<TH1D>("t0err", "t0err", 100, 0., 5.);
    t0_ = tfdir.make<TH1D>("t0", "t0", 339, 0., 1695.);
    d0_ = tfdir.make<TH1D>("d0", "d0", 200, -500., +500.);
    dout_ = tfdir.make<TH1D>("dout", "dout", 200, 0., +1000.);
    tanDip_ = tfdir.make<TH1D>("tanDip", "tanDip", 100, 0., 2.);
    momentum_ = tfdir.make<TH1D>("momentum", "momentum", 120, 100., 106.);
    radius_ = tfdir.make<TH1D>("radius", "radius", 160, 0., 800.);
    wavelength_ = tfdir.make<TH1D>("wavelength", "wavelength", 200, 0., 4000.);
    costh_ = tfdir.make<TH1D>("costh", "costh", 200, -1., +1.);
  }

  void HistTrackSum::fill(const TrackSummary& trk, double weight) {
    const auto& h = trk.states().at(0).helix();

    nactive_->Fill(trk.nactive(), weight);
    fitcon_->Fill(trk.fitcon(), weight);
    momerr_->Fill(trk.states().at(0).momentumError(), weight);
    t0err_->Fill(trk.t0Err(), weight);
    d0_->Fill(h.d0(), weight);
    dout_->Fill(h.dOut(), weight);
    t0_->Fill(trk.t0(), weight);
    tanDip_->Fill(h.tanDip(), weight);
    momentum_->Fill(trk.states().at(0).momentum().mag(), weight);
    radius_->Fill(h.radius(), weight);
    wavelength_->Fill(h.wavelength(), weight);
    costh_->Fill(trk.states().at(0).costh(), weight);
  }

}
