#include "Offline/Validation/inc/ValTrkCrvMatch.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include <cmath>
#include <map>

int mu2e::ValTrkCrvMatch::declare(const art::TFileDirectory& tfs) {
  GeomHandle<CosmicRayShield> crs;
  GeomHandle<KinKalGeom> kkg;
  auto const& shields = crs->getCRSScintillatorShields();
  int ns = static_cast<int>(shields.size());
  _v.assign(ns, KinKal::VEC3());
  _w.assign(ns, KinKal::VEC3());
  std::map<std::string, int> nameToShield;
  for (int i = 0; i < ns; ++i) nameToShield[shields[i].getName()] = i;
  for (auto const& sec : kkg->CRV()->sectors()) {
    auto it = nameToShield.find(sec.sname_);
    if (it == nameToShield.end()) continue;
    _v[it->second] = sec.sector_->vDirection();
    _w[it->second] = sec.sector_->normal();
  }

  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _hNTrk = tfs.make<TH1D>("NTrk", "N tracks", 21, -0.5, 20.5);
  _hNCand = tfs.make<TH1D>("NCand", "N candidates per track", 51, -0.5, 50.5);
  _hNGeom = tfs.make<TH1D>("NGeom", "N geometric candidates per track", 21, -0.5, 20.5);
  _hNXing = tfs.make<TH1D>("NXing", "N CRV crossings per track", 11, -0.5, 10.5);
  _hDtRef = tfs.make<TH1D>("dtRef", "t_{trk}(ref) - t_{cluster}, all candidates;ns", 200, -1000., 1000.);
  _hDtRefBest = tfs.make<TH1D>("dtRefBest", "t_{trk}(ref) - t_{cluster}, best candidate;ns", 200, -400., 400.);
  _hDtRefToF = tfs.make<TH1D>("dtRefToF", "dtRef - flight time, best candidate;ns", 200, -400., 400.);
  _hChi2 = tfs.make<TH1D>("chi2ndof", "#chi^{2}/ndof, best in sector", 100, 0., 20.);
  _hNShare = tfs.make<TH1D>("NShare", "tracks sharing the cluster", 11, -0.5, 10.5);
  _hFlags = tfs.make<TH1D>("Flags", "candidate flag bits", 16, -0.5, 15.5);
  _hXingShield = tfs.make<TH1D>("XingShield", "CRV crossings per shield", ns, -0.5, ns - 0.5);
  _hXingGeomShield = tfs.make<TH1D>("XingGeomShield", "crossings with a geometric candidate in their shield", ns, -0.5, ns - 0.5);
  _hDtXing = tfs.make<TH2D>("dtXing", "t_{cluster} - t_{trk} at the plane, best in sector;shield;ns", ns, -0.5, ns - 0.5, 100, -100., 100.);
  _hDu = tfs.make<TH2D>("du", "cluster - track along the bars, best in sector;shield;mm", ns, -0.5, ns - 0.5, 100, -3000., 3000.);
  _hDv = tfs.make<TH2D>("dv", "cluster - track across the bars, best in sector;shield;mm", ns, -0.5, ns - 0.5, 100, -500., 500.);
  _hPredMatched = tfs.make<TH2D>("PredMatched", "best geometric candidate;reference crossing shield;cluster shield", ns, -0.5, ns - 0.5, ns, -0.5, ns - 0.5);
  _pDvSlope.clear();
  for (int i = 0; i < ns; ++i) {
    std::string n = shields[i].getName();
    _pDvSlope.push_back(tfs.make<TProfile>(("dvVsSlope_" + n).c_str(), ("dv vs across-bar slope, " + n + ";slope;dv [mm]").c_str(), 40, -2., 2.));
  }
  _hProbePred = tfs.make<TH1D>("ProbePred", "layer probe: predicted bar per layer", 6, -0.5, 5.5);
  _hProbeHit = tfs.make<TH1D>("ProbeHit", "layer probe: pulse found per layer", 6, -0.5, 5.5);
  _hProbeInCluster = tfs.make<TH1D>("ProbeInCluster", "layer probe: pulse in a cluster per layer", 6, -0.5, 5.5);
  _hProbeDBar = tfs.make<TH1D>("ProbeDBar", "layer probe: selected - predicted bar", 21, -10.5, 10.5);
  _hProbeDt = tfs.make<TH1D>("ProbeDt", "layer probe: pulse - predicted time;ns", 100, -60., 60.);
  _hProbePE = tfs.make<TH1D>("ProbePE", "layer probe: PEs", 100, 0., 200.);
  _hProbeDu = tfs.make<TH1D>("ProbeDu", "layer probe: two-end position - predicted;mm", 100, -1500., 1500.);
  return 0;
}

int mu2e::ValTrkCrvMatch::fill(const mu2e::TrkCrvMatchCollection& coll, art::Event const& event) {
  // increment this by 1 any time the definitions of the histograms or the
  // histogram contents change, and will not match previous versions
  _hVer->Fill(1.0);

  int ns = static_cast<int>(_v.size());
  auto okShield = [&](int s) { return s >= 0 && s < ns; };
  _hNTrk->Fill(coll.size());
  for (auto const& m : coll) {
    _hNCand->Fill(m.cands_.size());
    _hNXing->Fill(m.xings_.size());
    int ngeom = 0;
    for (auto const& c : m.cands_) {
      _hDtRef->Fill(c.dtRef_);
      _hNShare->Fill(c.nShare_);
      for (int b = 0; b < 16; ++b)
        if (c.flags_ & (1u << b)) _hFlags->Fill(b);
      if (!c.hasFlag(TrkCrvCandidate::geomCand)) continue;
      ++ngeom;
      if (!c.hasFlag(TrkCrvCandidate::bestInSector) || !okShield(c.sector_)) continue;
      if (c.ndof_ > 0) _hChi2->Fill(c.chi2_ / c.ndof_);
      _hDtXing->Fill(c.sector_, c.dtXing_);
      if (!c.hasFlag(TrkCrvCandidate::singleEnded)) _hDu->Fill(c.sector_, c.du_);
      _hDv->Fill(c.sector_, c.dv_);
      if (c.xing_ >= 0 && static_cast<size_t>(c.xing_) < m.xings_.size()) {
        auto const& d = m.xings_[c.xing_].dir_;
        KinKal::VEC3 dir(d.X(), d.Y(), d.Z());
        double dw = dir.Dot(_w[c.sector_]);
        if (std::fabs(dw) > 1e-3) _pDvSlope[c.sector_]->Fill(dir.Dot(_v[c.sector_]) / dw, c.dv_);
      }
    }
    _hNGeom->Fill(ngeom);
    int best = m.bestGeom_ >= 0 ? m.bestGeom_ : m.bestTime_;
    if (best >= 0) {
      _hDtRefBest->Fill(m.cands_[best].dtRef_);
      _hDtRefToF->Fill(m.cands_[best].dtRefToF_);
    }
    if (m.bestGeom_ >= 0) {
      auto const& c = m.cands_[m.bestGeom_];
      if (c.xing_ >= 0 && static_cast<size_t>(c.xing_) < m.xings_.size()) _hPredMatched->Fill(m.xings_[c.xing_].sector_, c.sector_);
    }
    for (auto const& x : m.xings_) {
      _hXingShield->Fill(x.sector_);
      if (x.best_ >= 0) _hXingGeomShield->Fill(x.sector_);
      for (size_t L = 0; L < x.layers_.size(); ++L) {
        auto const& p = x.layers_[L];
        if (p.bar_ < 0) continue;
        _hProbePred->Fill(L);
        if (p.dBar_ == INT8_MAX) continue;
        _hProbeHit->Fill(L);
        if (p.inCluster_) _hProbeInCluster->Fill(L);
        _hProbeDBar->Fill(p.dBar_);
        _hProbeDt->Fill(p.dt_);
        _hProbePE->Fill(p.pe_);
        if (!std::isnan(p.uTdiv_)) _hProbeDu->Fill(p.uTdiv_ - p.uPred_);
      }
    }
  }
  return 0;
}
