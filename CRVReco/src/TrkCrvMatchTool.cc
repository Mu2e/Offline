#include "Offline/CRVReco/inc/TrkCrvMatchTool.hh"
#include "Offline/CRVReco/inc/TrkCrvPlaneError.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <algorithm>
#include <climits>
#include <cstdlib>

namespace mu2e {

  TrkCrvMatchTool::TrkCrvMatchTool(Config const& conf) :
    dtRefCenter_(conf.dtRefCenter()), maxDv_(conf.maxDv()), maxDu_(conf.maxDu()),
    doGeometric_(conf.doGeometric()), requireSameSector_(conf.requireSameSector()),
    sigmaV_(conf.sigmaV()), sigmaU_(conf.sigmaU()), sigmaT_(conf.sigmaT()),
    sigmaTSingle_(conf.sigmaTSingleEnded()), trkErrScale_(conf.trackErrorScale()), posErrQuad_(conf.trackPosErrQuad()),
    minTimeErr_(conf.minTrackTimeErr()), timeOffset_(conf.timeOffset()), keepAll_(conf.keepAll()),
    layerProbe_(conf.layerProbe()), probeBarWindow_(conf.probeBarWindow()), fiberSpeed_(conf.fiberSignalSpeed()),
    pulseDelay_(conf.pulseDelay())
  {
    auto w = conf.dtRefWindow(); auto x = conf.dtXingWindow(); auto p = conf.probeDtWindow();
    if(w.size() != 2 || x.size() != 2 || p.size() != 2) throw cet::exception("CONFIG") << "TrkCrvMatchTool: windows need 2 values";
    if(sigmaV_.empty() || sigmaU_.empty() || sigmaT_.empty()) throw cet::exception("CONFIG") << "TrkCrvMatchTool: empty resolution table";
    if(probeBarWindow_ < 0 || probeBarWindow_ > 100 || fiberSpeed_ <= 0) throw cet::exception("CONFIG") << "TrkCrvMatchTool: bad layer-probe settings";
    dtRefLo_ = w[0]; dtRefHi_ = w[1]; dtXingLo_ = x[0]; dtXingHi_ = x[1]; probeDtLo_ = p[0]; probeDtHi_ = p[1];
    for(auto const& s : conf.refSurfaces()) refSurfaces_.emplace_back(s);
  }

  double TrkCrvMatchTool::resolution(std::vector<double> const& tab, int sectorType) const {
    size_t i = sectorType > 0 ? static_cast<size_t>(sectorType - 1) : 0;
    return tab[std::min(i, tab.size() - 1)];
  }

  void TrkCrvMatchTool::beginRun() {
    GeomHandle<CosmicRayShield> crs;
    GeomHandle<DetectorSystem> det;
    GeomHandle<KinKalGeom> kkg;
    crs_ = crs.get();
    det_ = det.get();
    auto const& shields = crs_->getCRSScintillatorShields();
    sidToShield_.clear();
    sectors_.assign(shields.size(), SectorGeom());
    std::map<std::string,int> nameToShield;
    for(size_t i = 0; i < shields.size(); ++i) nameToShield[shields[i].getName()] = static_cast<int>(i);
    for(auto const& sec : kkg->CRV()->sectors()) {
      auto it = nameToShield.find(sec.sname_);
      if(it == nameToShield.end()) continue;
      SectorGeom& g = sectors_[it->second];
      g.valid = true;
      g.center = sec.sector_->center();
      g.u = sec.sector_->uDirection();
      g.v = sec.sector_->vDirection();
      g.w = sec.sector_->normal();
      g.uHalf = sec.sector_->uHalfLength();
      g.vHalf = sec.sector_->vHalfLength();
      try { sidToShield_[SurfaceId(sec.sname_)] = it->second; }
      catch(std::exception const&) { mf::LogWarning("TrkCrvMatchTool") << "no SurfaceId for CRV sector " << sec.sname_; }
    }
  }

  void TrkCrvMatchTool::setEvent(std::vector<art::Ptr<CrvCoincidenceCluster>> const& clusters, CrvRecoPulseCollection const* pulses) {
    if(!crs_) throw cet::exception("CONFIG") << "TrkCrvMatchTool: beginRun was not called";
    clusters_ = clusters;
    cinfo_.assign(clusters_.size(), ClusterInfo());
    clusteredPulses_.clear();
    for(size_t ic = 0; ic < clusters_.size(); ++ic) {
      auto const& cc = *clusters_[ic];
      auto& ci = cinfo_[ic];
      auto p = det_->toDetector(cc.GetAvgHitPos());
      ci.pos = KinKal::VEC3(p.x(), p.y(), p.z());
      ci.time = cc.GetAvgHitTime() + timeOffset_;
      ci.sectorType = cc.GetCrvSectorType();
      ci.singleEnded = !cc.HitPosAndTimeCalculated();
      std::map<int,double> pes;
      for(auto const& rp : cc.GetCrvRecoPulses()) {
        pes[crs_->getBar(rp->GetScintillatorBarIndex()).id().getShieldNumber()] += rp->GetPEs();
        clusteredPulses_.emplace(rp->GetScintillatorBarIndex().asInt(), rp->GetSiPMNumber(), rp->GetPulseTime());
      }
      double best = -1;
      for(auto const& [s, pe] : pes) if(pe > best) { best = pe; ci.shield = s; }
      ci.multiSector = pes.size() > 1;
    }
    pulses_ = layerProbe_ ? pulses : nullptr;
    pulsesByBar_.clear();
    if(pulses_) for(size_t ip = 0; ip < pulses_->size(); ++ip) pulsesByBar_[(*pulses_)[ip].GetScintillatorBarIndex().asInt()].push_back(ip);
  }

  bool TrkCrvMatchTool::sameCluster(art::Ptr<CrvCoincidenceCluster> const& a, art::Ptr<CrvCoincidenceCluster> const& b) {
    if(a.id() == b.id()) return a.key() == b.key();
    // a copy in another collection (e.g. SelectReco): same time, position and PEs
    return a->GetAvgHitTime() == b->GetAvgHitTime() && a->GetPEs() == b->GetPEs() && a->GetAvgHitPos() == b->GetAvgHitPos();
  }

  void TrkCrvMatchTool::probe(TrkCrvXing& x, double speed) const {
    if(x.sector_ < 0 || static_cast<size_t>(x.sector_) >= crs_->getCRSScintillatorShields().size()) return;
    auto const& shield = crs_->getCRSScintillatorShields()[x.sector_];
    auto const& bd = shield.getCRSScintillatorBarDetail();
    int td = bd.getThicknessDirection(), wd = bd.getWidthDirection(), ld = bd.getLengthDirection();
    double pitch = 2.0 * bd.getHalfWidth();
    CLHEP::Hep3Vector p0 = det_->toMu2e(CLHEP::Hep3Vector(x.pos_.X(), x.pos_.Y(), x.pos_.Z()));
    CLHEP::Hep3Vector d(x.dir_.X(), x.dir_.Y(), x.dir_.Z());
    auto const& modules = shield.getCRSScintillatorModules();
    if(modules.empty() || std::fabs(d[td]) < 1e-6) return;   // no layers, or a line parallel to them
    size_t nLayers = modules.front().getLayers().size();
    x.layers_.assign(nLayers, TrkCrvLayerProbe());
    for(size_t L = 0; L < nLayers; ++L) {
      auto& pr = x.layers_[L];
      // the predicted bar: nearest bar of the module whose layer L contains the line's crossing of that layer
      CRSScintillatorBar const* pbar = nullptr;
      CLHEP::Hep3Vector pp;
      double s = 0;
      for(auto const& mod : modules) {
        if(L >= mod.getLayers().size()) continue;
        auto const& lay = mod.getLayers()[L];
        double sl = (lay.getPosition()[td] - p0[td]) / d[td];
        CLHEP::Hep3Vector pl = p0 + sl * d;
        if(std::fabs(pl[wd] - lay.getPosition()[wd]) > lay.getHalfWidth() || std::fabs(pl[ld] - lay.getPosition()[ld]) > lay.getHalfLength()) continue;
        double bw = 1e30;
        for(auto const& b : lay.getBars()) {
          double dw = std::fabs(pl[wd] - b->getPosition()[wd]);
          if(dw < bw) { bw = dw; pbar = b.get(); }
        }
        pp = pl; s = sl;
        break;
      }
      if(!pbar) continue;
      pr.bar_ = static_cast<int16_t>(pbar->index().asInt());
      pr.uPred_ = pp[ld] - pbar->getPosition()[ld];
      pr.tPred_ = x.time_ + (speed > 0 ? s / speed : 0.0);
      if(!pulses_) continue;
      auto pulseDt = [&](CrvRecoPulse const& rp, CRSScintillatorBar const& b) {   // measured - predicted pulse time
        double dist = std::fabs(b.getSiPMPosition(rp.GetSiPMNumber())[ld] - pp[ld]);
        return rp.GetPulseTime() - pulseDelay_ + timeOffset_ - (pr.tPred_ + dist / fiberSpeed_);
      };
      // the selected bar: nearest to the predicted one within the bar window with an in-window pulse, then smallest |dt|
      CRSScintillatorBar const* sbar = nullptr;
      int sdb = 0, babs = INT_MAX;
      double bdt = 1e30;
      for(auto const& mod : modules) {
        if(L >= mod.getLayers().size()) continue;
        for(auto const& b : mod.getLayers()[L].getBars()) {
          int db = static_cast<int>(std::lround((b->getPosition()[wd] - pbar->getPosition()[wd]) / pitch));
          if(std::abs(db) > probeBarWindow_) continue;
          auto it = pulsesByBar_.find(b->index().asInt());
          if(it == pulsesByBar_.end()) continue;
          for(size_t ip : it->second) {
            double dt = pulseDt((*pulses_)[ip], *b);
            if(dt < probeDtLo_ || dt > probeDtHi_) continue;
            if(std::abs(db) < babs || (std::abs(db) == babs && std::fabs(dt) < std::fabs(bdt))) { babs = std::abs(db); bdt = dt; sbar = b.get(); sdb = db; }
          }
        }
      }
      if(!sbar) continue;
      pr.dBar_ = static_cast<int8_t>(sdb);
      pr.dt_ = bdt;
      double tside[2] = {NAN, NAN};
      for(size_t ip : pulsesByBar_.at(sbar->index().asInt())) {
        auto const& rp = (*pulses_)[ip];
        double dt = pulseDt(rp, *sbar);
        if(dt < probeDtLo_ || dt > probeDtHi_) continue;
        pr.pe_ += rp.GetPEs();
        if(clusteredPulses_.count({rp.GetScintillatorBarIndex().asInt(), rp.GetSiPMNumber(), rp.GetPulseTime()})) pr.inCluster_ = true;
        int side = sbar->getSiPMPosition(rp.GetSiPMNumber())[ld] < sbar->getPosition()[ld] ? 0 : 1;
        double t = rp.GetPulseTime();
        if(std::isnan(tside[side]) || t < tside[side]) tside[side] = t;
      }
      if(!std::isnan(tside[0]) && !std::isnan(tside[1])) pr.uTdiv_ = 0.5 * (tside[0] - tside[1]) * fiberSpeed_;
    }
  }

  TrkCrvMatch TrkCrvMatchTool::match(art::Ptr<KalSeed> const& seed, std::vector<art::Ptr<CrvCoincidenceCluster>> const& usedInFit) const {
    auto const& ks = *seed;
    TrkCrvMatch m;
    m.track_ = seed;
    double t0 = ks.t0Val();
    double t0err = std::sqrt(std::max(ks.t0Var(), 0.0));
    // reference times: the intersection on each surface nearest in time to t0
    auto nearest = [&](SurfaceId const& sid) -> int {
      int best = -1; double bdt = 1e30;
      for(size_t k = 0; k < ks.intersections().size(); ++k) {
        auto const& ki = ks.intersections()[k];
        if(ki.surfaceId() != sid) continue;
        double d = std::fabs(ki.time() - t0);
        if(d < bdt) { bdt = d; best = static_cast<int>(k); }
      }
      return best;
    };
    int kf = nearest(SurfaceId(SurfaceIdDetail::TT_Front));
    int km = nearest(SurfaceId(SurfaceIdDetail::TT_Mid));
    if(kf >= 0) { m.hasFront_ = true; m.tFront_ = ks.intersections()[kf].time(); }
    if(km >= 0) { m.hasMid_ = true; m.tMid_ = ks.intersections()[km].time(); }
    m.tRef_ = t0;
    bool hasRefPos = false;
    KinKal::VEC3 refPos;
    for(auto const& rs : refSurfaces_) {
      int k = nearest(rs);
      if(k < 0) continue;
      m.refSurface_ = rs;
      m.tRef_ = ks.intersections()[k].time();
      refPos = ks.intersections()[k].position3();
      hasRefPos = true;
      break;
    }
    m.timeOffset_ = timeOffset_;

    // CRV crossings
    double sigT = std::max(t0err, minTimeErr_);
    if(doGeometric_) {
      for(size_t k = 0; k < ks.intersections().size(); ++k) {
        auto const& ki = ks.intersections()[k];
        auto it = sidToShield_.find(ki.surfaceId());
        if(it == sidToShield_.end()) continue;
        auto const& g = sectors_[it->second];
        TrkCrvXing x;
        x.sid_ = ki.surfaceId();
        x.inter_ = static_cast<int16_t>(k);
        x.sector_ = static_cast<int16_t>(it->second);
        auto p = ki.position3();
        auto d = ki.momentum3().Unit();
        x.pos_ = XYZVectorF(p.X(), p.Y(), p.Z());
        x.dir_ = XYZVectorF(d.X(), d.Y(), d.Z());
        x.time_ = ki.time();
        x.inBounds_ = ki.inBounds();
        auto const& cov = ki.pstate_.stateCovariance();
        auto mom = ki.momentum3();
        x.sigU_ = trkPosErr(trkCrvPlaneError(cov, mom, g.u, g.w, 0.0));
        x.sigV_ = trkPosErr(trkCrvPlaneError(cov, mom, g.v, g.w, 0.0));
        x.sigT_ = sigT;
        x.mom_ = ki.mom();
        x.momErr_ = ki.momerr();
        KinKal::VEC3 loc = p - g.center;
        x.dEdgeU_ = g.uHalf - std::fabs(loc.Dot(g.u));
        x.dEdgeV_ = g.vHalf - std::fabs(loc.Dot(g.v));
        x.edgeSig_ = (x.sigU_ > 0 && x.sigV_ > 0) ? std::min(x.dEdgeU_ / x.sigU_, x.dEdgeV_ / x.sigV_) : 0.f;
        if(layerProbe_) probe(x, std::sqrt(ki.velocity().Mag2()));
        m.xings_.push_back(x);
      }
    }

    // candidates
    for(size_t ic = 0; ic < clusters_.size(); ++ic) {
      auto const& ci = cinfo_[ic];
      TrkCrvCandidate c;
      c.cluster_ = clusters_[ic];
      c.sector_ = static_cast<int16_t>(ci.shield);
      c.dtRef_ = m.tRef_ - ci.time;
      c.dtRefToF_ = hasRefPos ? c.dtRef_ - std::sqrt((refPos - ci.pos).Mag2()) / cLight_ : c.dtRef_;
      if(ci.singleEnded) c.setFlag(TrkCrvCandidate::singleEnded);
      if(ci.multiSector) c.setFlag(TrkCrvCandidate::multiSector);
      bool timeOK = keepAll_ || (c.dtRef_ >= dtRefLo_ && c.dtRef_ <= dtRefHi_);
      if(timeOK) c.setFlag(TrkCrvCandidate::timeCand);
      // geometric evidence: project from the CRV crossing nearest in path to the cluster's plane (straight line outside the field)
      int bestx = -1; double bds = 1e30;
      bool hasPlane = ci.shield >= 0 && static_cast<size_t>(ci.shield) < sectors_.size() && sectors_[ci.shield].valid;
      if(hasPlane) {
        auto const& g = sectors_[ci.shield];
        for(size_t k = 0; k < m.xings_.size(); ++k) {
          if(requireSameSector_ && m.xings_[k].sector_ != ci.shield) continue;
          auto const& xk = m.xings_[k];
          KinKal::VEC3 pk(xk.pos_.X(), xk.pos_.Y(), xk.pos_.Z()), dk(xk.dir_.X(), xk.dir_.Y(), xk.dir_.Z());
          double dn = dk.Dot(g.w);
          if(std::fabs(dn) < 1e-3) continue;
          double ds = std::fabs((ci.pos - pk).Dot(g.w) / dn);
          if(ds < bds) { bds = ds; bestx = static_cast<int>(k); }
        }
      }
      if(bestx >= 0) {
        auto const& g = sectors_[ci.shield];
        auto const& x = m.xings_[bestx];
        KinKal::VEC3 px(x.pos_.X(), x.pos_.Y(), x.pos_.Z()), dx(x.dir_.X(), x.dir_.Y(), x.dir_.Z());
        double dn = dx.Dot(g.w);
        auto const& ki = ks.intersections()[x.inter_];
        double speed = std::sqrt(ki.velocity().Mag2());
        double ds = (ci.pos - px).Dot(g.w) / dn;   // path length from the crossing to the cluster plane
        KinKal::VEC3 pp = px + ds * dx;
        double tp = x.time_ + (speed > 0 ? ds / speed : 0.0);
        KinKal::VEC3 r = ci.pos - pp;
        c.du_ = r.Dot(g.u);
        c.dv_ = r.Dot(g.v);
        c.dtXing_ = ci.time - tp;
        // track errors on the cluster's plane: the reference crossing's covariance carried over the straight step
        auto const& cov = ki.pstate_.stateCovariance();
        double tsu = trkPosErr(trkCrvPlaneError(cov, ki.momentum3(), g.u, g.w, ds));
        double tsv = trkPosErr(trkCrvPlaneError(cov, ki.momentum3(), g.v, g.w, ds));
        c.sigU_ = static_cast<float>(tsu);
        c.sigV_ = static_cast<float>(tsv);
        double sv = resolution(sigmaV_, ci.sectorType), su = resolution(sigmaU_, ci.sectorType);
        double st = ci.singleEnded ? sigmaTSingle_ : resolution(sigmaT_, ci.sectorType);
        c.chi2_ = c.dv_*c.dv_/(sv*sv + tsv*tsv) + c.dtXing_*c.dtXing_/(st*st + x.sigT_*x.sigT_);
        c.ndof_ = 2;
        if(!ci.singleEnded) { c.chi2_ += c.du_*c.du_/(su*su + tsu*tsu); c.ndof_ = 3; }
        bool geomOK = keepAll_ || (std::fabs(c.dv_) < maxDv_ && (ci.singleEnded || std::fabs(c.du_) < maxDu_)
                                   && c.dtXing_ >= dtXingLo_ && c.dtXing_ <= dtXingHi_);
        if(geomOK) {
          c.xing_ = static_cast<int16_t>(bestx);
          c.setFlag(TrkCrvCandidate::geomCand);
          if(x.sector_ == ci.shield) c.setFlag(TrkCrvCandidate::sameSector);
          KinKal::VEC3 loc = pp - g.center;
          if(std::fabs(loc.Dot(g.u)) > g.uHalf || std::fabs(loc.Dot(g.v)) > g.vHalf) c.setFlag(TrkCrvCandidate::outOfBounds);
        }
      }
      if(!c.hasFlag(TrkCrvCandidate::geomCand)) { c.chi2_ = -1; c.ndof_ = 0; c.du_ = c.dv_ = c.dtXing_ = c.sigU_ = c.sigV_ = 0; }
      if(c.hasFlag(TrkCrvCandidate::geomCand) || timeOK) {
        for(auto const& u : usedInFit) if(sameCluster(u, c.cluster_)) { c.setFlag(TrkCrvCandidate::usedInFit); break; }
        m.cands_.push_back(c);
      }
    }

    // rank: geometric by chi2, then time-only by distance from the time centre
    std::stable_sort(m.cands_.begin(), m.cands_.end(), [&](TrkCrvCandidate const& a, TrkCrvCandidate const& b) {
      bool ga = a.hasFlag(TrkCrvCandidate::geomCand), gb = b.hasFlag(TrkCrvCandidate::geomCand);
      if(ga != gb) return ga;
      if(ga) return a.chi2_ < b.chi2_;
      return std::fabs(a.dtRef_ - dtRefCenter_) < std::fabs(b.dtRef_ - dtRefCenter_);
    });
    double bt = 1e30;
    std::map<int,int16_t> bestInSector;   // cluster sector -> lowest-chi2 geometric candidate
    for(size_t j = 0; j < m.cands_.size(); ++j) {
      auto& c = m.cands_[j];
      if(c.hasFlag(TrkCrvCandidate::geomCand)) {
        if(m.bestGeom_ < 0) m.bestGeom_ = static_cast<int16_t>(j);
        if(bestInSector.emplace(c.sector_, static_cast<int16_t>(j)).second) c.setFlag(TrkCrvCandidate::bestInSector);
      }
      if(c.hasFlag(TrkCrvCandidate::timeCand) && std::fabs(c.dtRef_ - dtRefCenter_) < bt) {
        bt = std::fabs(c.dtRef_ - dtRefCenter_);
        m.bestTime_ = static_cast<int16_t>(j);
      }
    }
    for(auto& x : m.xings_) {
      auto ib = bestInSector.find(x.sector_);
      if(ib != bestInSector.end()) x.best_ = ib->second;
    }
    int bf = m.bestGeom_ >= 0 ? m.bestGeom_ : m.bestTime_;
    if(bf >= 0) m.cands_[bf].setFlag(TrkCrvCandidate::bestForTrack);
    return m;
  }

  void TrkCrvMatchTool::resolve(TrkCrvMatchCollection& out) const {
    struct Ref { size_t trk, cand; };
    std::map<size_t, std::vector<Ref>> byCluster;
    for(size_t it = 0; it < out.size(); ++it)
      for(size_t j = 0; j < out[it].cands_.size(); ++j) byCluster[out[it].cands_[j].cluster_.key()].push_back({it, j});
    auto metric = [&](Ref const& r) {   // geometric chi2 first, then time distance (offset so geometric always wins)
      auto const& c = out[r.trk].cands_[r.cand];
      return c.hasFlag(TrkCrvCandidate::geomCand) ? static_cast<double>(c.chi2_) : 1e9 + std::fabs(c.dtRef_ - dtRefCenter_);
    };
    std::vector<std::pair<double,Ref>> geomPairs;
    for(auto& [key, refs] : byCluster) {
      Ref best = refs.front();
      for(auto const& r : refs) {
        out[r.trk].cands_[r.cand].nShare_ = static_cast<uint8_t>(std::min<size_t>(refs.size(), 255));
        if(metric(r) < metric(best)) best = r;
        if(out[r.trk].cands_[r.cand].hasFlag(TrkCrvCandidate::geomCand)) geomPairs.push_back({metric(r), r});
      }
      out[best.trk].cands_[best.cand].setFlag(TrkCrvCandidate::bestForCluster);
    }
    std::stable_sort(geomPairs.begin(), geomPairs.end(), [](auto const& a, auto const& b) { return a.first < b.first; });
    std::map<size_t,bool> usedCluster;
    std::map<std::pair<size_t,int>,bool> usedXing;   // (track, cluster sector)
    for(auto const& [chi2, r] : geomPairs) {
      auto& c = out[r.trk].cands_[r.cand];
      auto xk = std::make_pair(r.trk, static_cast<int>(c.sector_));
      if(usedCluster[c.cluster_.key()] || usedXing[xk]) continue;
      usedCluster[c.cluster_.key()] = usedXing[xk] = true;
      c.setFlag(TrkCrvCandidate::exclusive);
    }
  }
}
