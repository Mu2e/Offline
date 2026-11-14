#include "Offline/Validation/inc/ValTrkCrvMatch.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include <cmath>

int mu2e::ValTrkCrvMatch::declare(const art::TFileDirectory& tfs) {
  int ns = _axis.n();
  double lo = _axis.lo(), hi = _axis.hi();
  _hVer = tfs.make<TH1D>("Ver", "Version Number", 101, -0.5, 100.0);
  _hNTrk = tfs.make<TH1D>("NTrk", "N tracks", 21, -0.5, 20.5);
  _hNCand = tfs.make<TH1D>("NCand", "N candidates per track", 51, -0.5, 50.5);
  _hNGeom = tfs.make<TH1D>("NGeom", "N geometric candidates per track", 21, -0.5, 20.5);
  _hNXing = tfs.make<TH1D>("NXing", "N CRV crossings per track (extrapolated, in bounds)", 11, -0.5, 10.5);
  _hNXingProj = tfs.make<TH1D>("NXingProj", "N projected CRV crossings per line track with no crossing", 11, -0.5, 10.5);
  _hDtRef = tfs.make<TH1D>("dtRef", "t_{trk}(ref) - t_{cluster}, all candidates;ns", 200, -1000., 1000.);
  _hDtRefBest = tfs.make<TH1D>("dtRefBest", "t_{trk}(ref) - t_{cluster}, best candidate;ns", 200, -400., 400.);
  _hDtRefToF = tfs.make<TH1D>("dtRefToF", "dtRef - flight time, best candidate;ns", 200, -400., 400.);
  _hChi2 = tfs.make<TH1D>("chi2ndof", "#chi^{2}/ndof, best in sector, from an in-bounds crossing", 100, 0., 20.);
  _hChi2Proj = tfs.make<TH1D>("chi2ndofProj", "#chi^{2}/ndof, best in sector, from a projected crossing", 100, 0., 20.);
  _hNShare = tfs.make<TH1D>("NShare", "tracks sharing the cluster", 11, -0.5, 10.5);
  _hFlags = tfs.make<TH1D>("Flags", "candidate flag bits", 16, -0.5, 15.5);
  _hXingShield = tfs.make<TH1D>("XingShield", "CRV crossings per sector (extrapolated, in bounds)", ns, lo, hi);
  _hXingGeomShield = tfs.make<TH1D>("XingGeomShield", "in-bounds crossings with a geometric candidate in their sector", ns, lo, hi);
  _hXingProjShield = tfs.make<TH1D>("XingProjShield", "projected crossings per sector (line tracks with no crossing)", ns, lo, hi);
  _hDtXing = tfs.make<TH2D>("dtXing", "t_{cluster} - t_{trk} at the plane, best in sector;cluster sector;ns", ns, lo, hi, 100, -100., 100.);
  _hDu = tfs.make<TH2D>("du", "cluster - track along the bars, best in sector;cluster sector;mm", ns, lo, hi, 100, -3000., 3000.);
  _hDv = tfs.make<TH2D>("dv", "cluster - track across the bars, best in sector;cluster sector;mm", ns, lo, hi, 100, -500., 500.);
  _hPredMatched = tfs.make<TH2D>("PredMatched", "best geometric candidate;reference crossing sector;cluster sector", ns, lo, hi, ns, lo, hi);
  _hSlopeNoGeom = tfs.make<TH1D>("SlopeNoGeom", "best in sector, not in dvVsSlope: sector not in this job's geometry;cluster sector", ns, lo, hi);
  for (TH1* h : {static_cast<TH1*>(_hXingShield), static_cast<TH1*>(_hXingGeomShield), static_cast<TH1*>(_hXingProjShield), static_cast<TH1*>(_hDtXing),
                 static_cast<TH1*>(_hDu), static_cast<TH1*>(_hDv), static_cast<TH1*>(_hPredMatched), static_cast<TH1*>(_hSlopeNoGeom)})
    _axis.label(h->GetXaxis());
  _axis.label(_hPredMatched->GetYaxis());

  // dv vs slope needs the sector directions, from this job's geometry: one profile per sector it has
  GeomHandle<KinKalGeom> kkg;
  for (auto const& sec : kkg->CRV()->sectors()) {
    int b = _axis.bin(sec.sname_);
    if (b == _axis.unknown() || _plane.count(b)) continue;
    auto& pl = _plane[b];
    pl.c = sec.sector_->center();
    pl.v = sec.sector_->vDirection();
    pl.w = sec.sector_->normal();
    pl.p = tfs.make<TProfile>(("dvVsSlope_" + sec.sname_).c_str(), ("dv vs across-bar slope, " + sec.sname_ + ";slope;dv [mm]").c_str(), 40, -2., 2.);
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
  _hVer->Fill(2.0);

  // a cluster's sector: the SurfaceId of a crossing with the same shield index in the product's geometry
  auto clusterBin = [&](int16_t shield) {
    auto it = _shieldBin.find(shield);
    return it == _shieldBin.end() ? _axis.unknown() : it->second;
  };
  _hNTrk->Fill(coll.size());
  for (auto const& m : coll) {
    for (auto const& x : m.xings_) _shieldBin[x.sector_] = _axis.bin(x.sid_);
    _hNCand->Fill(m.cands_.size());
    int nproj = 0;
    for (auto const& x : m.xings_) nproj += !x.inBounds_;
    _hNXing->Fill(m.xings_.size() - nproj);
    if (nproj > 0) _hNXingProj->Fill(nproj);
    int ngeom = 0;
    for (auto const& c : m.cands_) {
      _hDtRef->Fill(c.dtRef_);
      _hNShare->Fill(c.nShare_);
      for (int b = 0; b < 16; ++b)
        if (c.flags_ & (1u << b)) _hFlags->Fill(b);
      if (!c.hasFlag(TrkCrvCandidate::geomCand)) continue;
      ++ngeom;
      if (!c.hasFlag(TrkCrvCandidate::bestInSector)) continue;
      int b = clusterBin(c.sector_);
      bool proj = c.xing_ >= 0 && static_cast<size_t>(c.xing_) < m.xings_.size() && !m.xings_[c.xing_].inBounds_;
      if (c.ndof_ > 0) (proj ? _hChi2Proj : _hChi2)->Fill(c.chi2_ / c.ndof_);
      _hDtXing->Fill(b, c.dtXing_);
      if (!c.hasFlag(TrkCrvCandidate::singleEnded)) _hDu->Fill(b, c.du_);
      _hDv->Fill(b, c.dv_);
      if (c.xing_ < 0 || static_cast<size_t>(c.xing_) >= m.xings_.size()) continue;
      // the job's directions apply only if its plane of the reference crossing's sector holds that crossing
      auto const& x = m.xings_[c.xing_];
      auto ix = _plane.find(_axis.bin(x.sid_));
      auto ic = _plane.find(b);
      KinKal::VEC3 pos(x.pos_.X(), x.pos_.Y(), x.pos_.Z());
      if (ix == _plane.end() || ic == _plane.end() || std::fabs((pos - ix->second.c).Dot(ix->second.w)) > _planeTol) {
        _hSlopeNoGeom->Fill(b);
        continue;
      }
      KinKal::VEC3 dir(x.dir_.X(), x.dir_.Y(), x.dir_.Z());
      double dw = dir.Dot(ic->second.w);
      if (std::fabs(dw) > 1e-3) ic->second.p->Fill(dir.Dot(ic->second.v) / dw, c.dv_);
    }
    _hNGeom->Fill(ngeom);
    int best = m.bestGeom_ >= 0 ? m.bestGeom_ : m.bestTime_;
    if (best >= 0) {
      _hDtRefBest->Fill(m.cands_[best].dtRef_);
      _hDtRefToF->Fill(m.cands_[best].dtRefToF_);
    }
    if (m.bestGeom_ >= 0) {
      auto const& c = m.cands_[m.bestGeom_];
      if (c.xing_ >= 0 && static_cast<size_t>(c.xing_) < m.xings_.size())
        _hPredMatched->Fill(_axis.bin(m.xings_[c.xing_].sid_), clusterBin(c.sector_));
    }
    for (auto const& x : m.xings_) {
      if (!x.inBounds_) {
        _hXingProjShield->Fill(_axis.bin(x.sid_));
        continue;
      }
      _hXingShield->Fill(_axis.bin(x.sid_));
      if (x.best_ >= 0) _hXingGeomShield->Fill(_axis.bin(x.sid_));
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
