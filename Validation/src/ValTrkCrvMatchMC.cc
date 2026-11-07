#include "Offline/Validation/inc/ValTrkCrvMatchMC.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include <cmath>

int mu2e::ValTrkCrvMatchMC::declare(const art::TFileDirectory& tfs) {
  GeomHandle<CosmicRayShield> crs;
  _ns = static_cast<int>(crs->getCRSScintillatorShields().size());
  double lo = -0.5, hi = _ns - 0.5;
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _hBestGeomRel = tfs.make<TH1D>("BestGeomRel", "best geometric candidate, MCRelationship to the track's particle: 0 no MC, 1 unrelated, 2 same, 3 daughter, 4 mother, 5 sibling, 6 udaughter, 7 umother, 8 usibling", 9, -0.5, 8.5);
  _hBestTimeRel = tfs.make<TH1D>("BestTimeRel", "best time candidate, MCRelationship to the track's particle", 9, -0.5, 8.5);
  _hGeomCandRel = tfs.make<TH1D>("GeomCandRel", "all geometric candidates, MCRelationship", 9, -0.5, 8.5);
  _hTimeCandRel = tfs.make<TH1D>("TimeCandRel", "time-only candidates, MCRelationship", 9, -0.5, 8.5);
  _hTrueClu = tfs.make<TH1D>("TrueClu", "own cluster of the track's particle per shield", _ns, lo, hi);
  _hTrueCluGeom = tfs.make<TH1D>("TrueCluGeom", "own cluster is a geometric candidate", _ns, lo, hi);
  _hTrueCluBest = tfs.make<TH1D>("TrueCluBest", "own cluster is the best in its sector", _ns, lo, hi);
  _hRecoXing = tfs.make<TH1D>("RecoXing", "reco crossings per shield", _ns, lo, hi);
  _hRecoXingTrue = tfs.make<TH1D>("RecoXingTrue", "reco crossings with a true crossing of the same shield", _ns, lo, hi);
  _hTrueXing = tfs.make<TH1D>("TrueXing", "true crossings of the track's particle per shield", _ns, lo, hi);
  _hTrueXingFound = tfs.make<TH1D>("TrueXingFound", "true crossings with a reco crossing", _ns, lo, hi);
  _hXingDist = tfs.make<TH2D>("XingDist", "|reco - true crossing|;shield;mm", _ns, lo, hi, 100, 0., 500.);
  _hXingDt = tfs.make<TH2D>("XingDt", "reco - true crossing time;shield;ns", _ns, lo, hi, 100, -20., 20.);
  return 0;
}

int mu2e::ValTrkCrvMatchMC::fill(const mu2e::TrkCrvMatchMCCollection& coll, art::Event const& event) {
  // increment this by 1 any time the definitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(1.0);

  auto ok = [&](int s) { return s >= 0 && s < _ns; };
  for (auto const& mc : coll) {
    if (mc.match_.isNull() || !mc.match_.isAvailable()) continue;
    auto const& m = *mc.match_;
    if (mc.cands_.size() != m.cands_.size()) continue;
    auto relBin = [](TrkCrvCandidateMC const& cm) { return cm.hasMC_ ? cm.rel_.relationship() + 2 : 0; };   // 0 no MC, 1 none, 2 same, ...
    if (m.bestGeom_ >= 0) _hBestGeomRel->Fill(relBin(mc.cands_[m.bestGeom_]));
    if (m.bestTime_ >= 0) _hBestTimeRel->Fill(relBin(mc.cands_[m.bestTime_]));
    for (size_t j = 0; j < m.cands_.size(); ++j)
      (m.cands_[j].hasFlag(TrkCrvCandidate::geomCand) ? _hGeomCandRel : _hTimeCandRel)->Fill(relBin(mc.cands_[j]));
    for (auto const& xm : mc.xings_) {
      if (xm.rel_ != MCRelationship::same || !ok(xm.sector_)) continue;
      if (xm.geom_) {
        _hTrueXing->Fill(xm.sector_);
        if (xm.recoXing_ >= 0) _hTrueXingFound->Fill(xm.sector_);
      }
      if (xm.trueCluster_ < 0) continue;
      _hTrueClu->Fill(xm.sector_);
      if (xm.trueCand_ < 0 || static_cast<size_t>(xm.trueCand_) >= m.cands_.size()) continue;
      auto const& c = m.cands_[xm.trueCand_];
      if (!c.hasFlag(TrkCrvCandidate::geomCand)) continue;
      _hTrueCluGeom->Fill(xm.sector_);
      if (c.hasFlag(TrkCrvCandidate::bestInSector)) _hTrueCluBest->Fill(xm.sector_);
    }
    for (size_t ir = 0; ir < m.xings_.size() && ir < mc.recoXingTruth_.size(); ++ir) {
      auto const& rx = m.xings_[ir];
      if (!ok(rx.sector_)) continue;
      _hRecoXing->Fill(rx.sector_);
      int t = mc.recoXingTruth_[ir];
      if (t < 0 || static_cast<size_t>(t) >= mc.xings_.size()) continue;
      _hRecoXingTrue->Fill(rx.sector_);
      auto const& xm = mc.xings_[t];
      if (!xm.geom_ || xm.rel_ != MCRelationship::same) continue;
      _hXingDist->Fill(rx.sector_, std::sqrt((rx.pos_ - xm.pos_).Mag2()));
      _hXingDt->Fill(rx.sector_, rx.time_ - xm.time_);
    }
  }
  return 0;
}
