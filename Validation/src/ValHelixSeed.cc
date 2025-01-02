#include "Offline/Validation/inc/ValHelixSeed.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/Mu2eUtilities/inc/HelixTool.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include <cmath>

int mu2e::ValHelixSeed::declare(const art::TFileDirectory& tfs) {
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.5);
  _hN = tfs.make<TH1D>("NSeed", "N KalSeed", 11, -0.5, 10.5);
  _hNCombo = tfs.make<TH1D>("NCombo", "N Combo Hits", 101, -0.5, 100.5);
  _hNStrHit = tfs.make<TH1D>("NStrHit", "N Straw Hits", 101, -0.5, 100.5);
  _hStatus = tfs.make<TH1D>("Status", "Status", 32, -0.5, 31.5);
  _ht0 = tfs.make<TH1D>("t0", "t0", 100, 400.0, 1800.0);
  _hp = tfs.make<TH1D>("p", "p (mm)", 200, -400.0, 400.);
  _hpce = tfs.make<TH1D>("pce", "p CE (mm)", 100, 300.0, 400.);
  _hpt = tfs.make<TH1D>("pt", "pt (mm)", 100, 0., 400.);
  _hD0 = tfs.make<TH1D>("d0", "d0", 400, -400., 400.);
  _hPhi0 = tfs.make<TH1D>("phi0", "phi0", 100, -M_PI, M_PI);
  _hLambda = tfs.make<TH1D>("lambda", "lambda", 100, -400.0, 400.0);
  _hchi2dXY = tfs.make<TH1D>("chi2dXY", "chi2dXY", 100, 0.0, 10.0);
  _hchi2dZPhi = tfs.make<TH1D>("chi2dZPhi", "chi2dZPhi", 100, 0.0, 10.0);

  return 0;
}

int mu2e::ValHelixSeed::fill(const mu2e::HelixSeedCollection& coll,
                             art::Event const& event) {
  // increment this by 1 any time the defnitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(0.0);
  _hN->Fill(coll.size());
  mu2e::GeomHandle<mu2e::Tracker> th;
  auto myTracker = th.get();

  for (auto const& hs : coll) {
    _hNCombo->Fill(hs.hits().size());
    HelixTool helTool(&hs, myTracker);
    int nstrawhits = helTool.nstrawhits();
    _hNStrHit->Fill(nstrawhits);
    const TrkFitFlag& tff = hs.status();

    for (auto sn : tff.bitNames()) {
      if (tff.hasAnyProperty(TrkFitFlag(sn.first)))
        _hStatus->Fill(std::log2(sn.second));
    }

    _ht0->Fill(hs.t0().t0());

    auto const& rh = hs.helix();

    double p = rh.momentum();
    double charge(1.);

    if ( ((rh.lambda()>0.) && (hs.recoDir().slope()>0.)) ||
         ((rh.lambda()<0.) && (hs.recoDir().slope()<0.)) ){
      charge = -1.;
    }

    _hp->Fill(p*charge);
    _hpce->Fill(p);
    _hpt->Fill(rh.radius());
    _hD0->Fill(rh.rcent() - rh.radius());
    _hPhi0->Fill(rh.fz0());
    _hLambda->Fill(rh.lambda());
    _hchi2dXY->Fill(rh.chi2dXY());
    _hchi2dZPhi->Fill(rh.chi2dZPhi());
  }
  return 0;
}
