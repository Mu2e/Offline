#include "Offline/Validation/inc/ValTrkCrvMatchMC.hh"
#include <cmath>

int mu2e::ValTrkCrvMatchMC::declare(const art::TFileDirectory& tfs) {
  int ns = _axis.n();
  double lo = _axis.lo(), hi = _axis.hi();
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _hBestGeomRel = tfs.make<TH1D>("BestGeomRel", "best geometric candidate, MCRelationship to the track's particle: 0 no MC, 1 unrelated, 2 same, 3 daughter, 4 mother, 5 sibling, 6 udaughter, 7 umother, 8 usibling", 9, -0.5, 8.5);
  _hBestTimeRel = tfs.make<TH1D>("BestTimeRel", "best time candidate, MCRelationship to the track's particle", 9, -0.5, 8.5);
  _hGeomCandRel = tfs.make<TH1D>("GeomCandRel", "all geometric candidates, MCRelationship", 9, -0.5, 8.5);
  _hTimeCandRel = tfs.make<TH1D>("TimeCandRel", "time-only candidates, MCRelationship", 9, -0.5, 8.5);
  _hTrueClu = tfs.make<TH1D>("TrueClu", "own cluster of the track's particle per sector", ns, lo, hi);
  _hTrueCluGeom = tfs.make<TH1D>("TrueCluGeom", "own cluster is a geometric candidate", ns, lo, hi);
  _hTrueCluBest = tfs.make<TH1D>("TrueCluBest", "own cluster is the best in its sector", ns, lo, hi);
  _hRecoXing = tfs.make<TH1D>("RecoXing", "reco crossings per sector (extrapolated, in bounds)", ns, lo, hi);
  _hRecoXingTrue = tfs.make<TH1D>("RecoXingTrue", "in-bounds reco crossings with a true crossing of the same sector", ns, lo, hi);
  _hRecoXingProj = tfs.make<TH1D>("RecoXingProj", "projected crossings per sector (line tracks with no crossing)", ns, lo, hi);
  _hRecoXingProjTrue = tfs.make<TH1D>("RecoXingProjTrue", "projected crossings with a true crossing of the same sector", ns, lo, hi);
  _hTrueXing = tfs.make<TH1D>("TrueXing", "true crossings of the track's particle per sector", ns, lo, hi);
  _hTrueXingFound = tfs.make<TH1D>("TrueXingFound", "true crossings with an in-bounds reco crossing", ns, lo, hi);
  _hTrueXingProj = tfs.make<TH1D>("TrueXingProj", "true crossings matched only by a projected crossing", ns, lo, hi);
  _hXingDist = tfs.make<TH2D>("XingDist", "|reco - true crossing|, in-bounds reco;sector;mm", ns, lo, hi, 100, 0., 500.);
  _hXingDistProj = tfs.make<TH2D>("XingDistProj", "|reco - true crossing|, projected reco;sector;mm", ns, lo, hi, 100, 0., 500.);
  _hXingDt = tfs.make<TH2D>("XingDt", "reco - true crossing time;sector;ns", ns, lo, hi, 100, -20., 20.);
  for (TH1* h : {static_cast<TH1*>(_hTrueClu), static_cast<TH1*>(_hTrueCluGeom), static_cast<TH1*>(_hTrueCluBest),
                 static_cast<TH1*>(_hRecoXing), static_cast<TH1*>(_hRecoXingTrue), static_cast<TH1*>(_hRecoXingProj),
                 static_cast<TH1*>(_hRecoXingProjTrue), static_cast<TH1*>(_hTrueXing), static_cast<TH1*>(_hTrueXingFound),
                 static_cast<TH1*>(_hTrueXingProj), static_cast<TH1*>(_hXingDist), static_cast<TH1*>(_hXingDistProj),
                 static_cast<TH1*>(_hXingDt)})
    _axis.label(h->GetXaxis());
  return 0;
}

int mu2e::ValTrkCrvMatchMC::fill(const mu2e::TrkCrvMatchMCCollection& coll, art::Event const& event) {
  // increment this by 1 any time the definitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(2.0);

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
      if (xm.rel_ != MCRelationship::same) continue;
      int b = _axis.bin(xm.sid_);
      if (xm.geom_) {
        _hTrueXing->Fill(b);
        bool reco = xm.recoXing_ >= 0 && static_cast<size_t>(xm.recoXing_) < m.xings_.size();
        if (reco) (m.xings_[xm.recoXing_].inBounds_ ? _hTrueXingFound : _hTrueXingProj)->Fill(b);
      }
      if (xm.trueCluster_ < 0) continue;
      _hTrueClu->Fill(b);
      if (xm.trueCand_ < 0 || static_cast<size_t>(xm.trueCand_) >= m.cands_.size()) continue;
      auto const& c = m.cands_[xm.trueCand_];
      if (!c.hasFlag(TrkCrvCandidate::geomCand)) continue;
      _hTrueCluGeom->Fill(b);
      if (c.hasFlag(TrkCrvCandidate::bestInSector)) _hTrueCluBest->Fill(b);
    }
    for (size_t ir = 0; ir < m.xings_.size() && ir < mc.recoXingTruth_.size(); ++ir) {
      auto const& rx = m.xings_[ir];
      int b = _axis.bin(rx.sid_);
      bool proj = !rx.inBounds_;
      (proj ? _hRecoXingProj : _hRecoXing)->Fill(b);
      int t = mc.recoXingTruth_[ir];
      if (t < 0 || static_cast<size_t>(t) >= mc.xings_.size()) continue;
      (proj ? _hRecoXingProjTrue : _hRecoXingTrue)->Fill(b);
      auto const& xm = mc.xings_[t];
      if (!xm.geom_ || xm.rel_ != MCRelationship::same) continue;
      (proj ? _hXingDistProj : _hXingDist)->Fill(b, std::sqrt((rx.pos_ - xm.pos_).Mag2()));
      if (!proj) _hXingDt->Fill(b, rx.time_ - xm.time_);
    }
  }
  return 0;
}
