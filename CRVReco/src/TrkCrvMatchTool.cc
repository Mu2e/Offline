#include "Offline/CRVReco/inc/TrkCrvMatchTool.hh"
#include "Offline/CRVReco/inc/TrkCrvPlaneError.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "cetlib_except/exception.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <algorithm>
#include <climits>
#include <cstdlib>

namespace mu2e {

  TrkCrvMatchTool::TrkCrvMatchTool(Config const& conf) :
    earliestRef_(conf.earliestReference()), dtRefCenter_(conf.dtRefCenter()), maxDv_(conf.maxDv()), maxDu_(conf.maxDu()),
    doGeometric_(conf.doGeometric()), requireSameSector_(conf.requireSameSector()),
    sigmaV_(conf.sigmaV()), sigmaU_(conf.sigmaU()), sigmaT_(conf.sigmaT()),
    sigmaTSingle_(conf.sigmaTSingleEnded()), trkErrScale_(conf.trackErrorScale()), posErrQuad_(conf.trackPosErrQuad()),
    minTimeErr_(conf.minTrackTimeErr()), timeOffset_(conf.timeOffset()), keepAll_(conf.keepAll()),
    layerProbe_(conf.layerProbe()), probeAllPulses_(conf.probeAllPulses()), probeBarWindow_(conf.probeBarWindow()), fiberSpeed_(conf.fiberSignalSpeed()),
    pulseDelay_(conf.pulseDelay())
  {
    auto w = conf.dtRefWindow(); auto x = conf.dtXingWindow(); auto p = conf.probeDtWindow();
    if(w.size() != 2 || x.size() != 2 || p.size() != 2) throw cet::exception("CONFIG") << "TrkCrvMatchTool: windows need 2 values";
    if(sigmaV_.empty() || sigmaU_.empty() || sigmaT_.empty()) throw cet::exception("CONFIG") << "TrkCrvMatchTool: empty resolution table";
    if(probeBarWindow_ < 0 || probeBarWindow_ > 100 || fiberSpeed_ <= 0) throw cet::exception("CONFIG") << "TrkCrvMatchTool: bad layer-probe settings";
    dtRefLo_ = w[0]; dtRefHi_ = w[1]; dtXingLo_ = x[0]; dtXingHi_ = x[1]; probeDtLo_ = p[0]; probeDtHi_ = p[1];
    for(auto const& s : conf.refSurfaces()) refSurfaces_.emplace_back(s);
  }

  TrkCrvPlaneCov TrkCrvMatchTool::errorModel(TrkCrvCov2 const& cov) const {
    TrkCrvPlaneCov out;
    for(int i = 0; i < 2; ++i)
      for(int j = i; j < 2; ++j) out(i, j) = trkErrScale_ * trkErrScale_ * cov(i, j) + (i == j ? posErrQuad_ * posErrQuad_ : 0.0);
    return out;
  }

  double TrkCrvMatchTool::resolution(std::vector<double> const& tab, int sectorType) const {
    size_t i = sectorType > 0 ? static_cast<size_t>(sectorType - 1) : 0;
    return tab[std::min(i, tab.size() - 1)];
  }

  void TrkCrvMatchTool::beginRun() {
    GeomHandle<CosmicRayShield> crs;
    GeomHandle<DetectorSystem> det;
    crs_ = crs.get();
    det_ = det.get();
    sectors_ = trkCrvSectors();
    sidToShield_.clear();
    for(size_t i = 0; i < sectors_.size(); ++i)
      if(sectors_[i].valid && sectors_[i].sid != SurfaceId()) sidToShield_[sectors_[i].sid] = static_cast<int>(i);
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
      ci.shield = trkCrvDominantShield(cc, *crs_, &ci.multiSector);
      for(auto const& rp : cc.GetCrvRecoPulses())
        clusteredPulses_.emplace(rp->GetScintillatorBarIndex().asInt(), rp->GetSiPMNumber(), rp->GetPulseTime());
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

  bool TrkCrvMatchTool::projectXing(KalSeed const& ks, int shield, double sigT, TrkCrvXing& x) const {
    auto const& g = sectors_[shield];
    // the straight line from the KalIntersection nearest in path to the sector's extended mid-plane
    int kbest = -1;
    double sbest = 0;
    for(size_t k = 0; k < ks.intersections().size(); ++k) {
      auto const& ki = ks.intersections()[k];
      double dn = ki.momentum3().Unit().Dot(g.w);
      if(std::fabs(dn) < 1e-3) continue;
      double s = (g.center - ki.position3()).Dot(g.w) / dn;
      if(kbest < 0 || std::fabs(s) < std::fabs(sbest)) { kbest = static_cast<int>(k); sbest = s; }
    }
    if(kbest < 0) return false;
    auto const& ki = ks.intersections()[kbest];
    auto d = ki.momentum3().Unit();
    KinKal::VEC3 p = ki.position3() + sbest * d;
    KinKal::VEC3 loc = p - g.center;
    double dEdgeU = g.uHalf - std::fabs(loc.Dot(g.u)), dEdgeV = g.vHalf - std::fabs(loc.Dot(g.v));
    if(dEdgeU >= 0 && dEdgeV >= 0) return false;   // in-bounds crossings belong to the extrapolation's intersections
    double speed = std::sqrt(ki.velocity().Mag2());
    x.sid_ = g.sid;
    x.inter_ = static_cast<int16_t>(kbest);
    x.sector_ = static_cast<int16_t>(shield);
    x.pos_ = XYZVectorF(p.X(), p.Y(), p.Z());
    x.dir_ = XYZVectorF(d.X(), d.Y(), d.Z());
    x.time_ = ki.time() + (speed > 0 ? sbest / speed : 0.0);
    x.inBounds_ = false;
    x.cov_ = errorModel(trkCrvPlaneCovariance(ki.pstate_.stateCovariance(), ki.momentum3(), g.u, g.v, g.w, sbest));
    x.sigT_ = sigT;
    x.mom_ = ki.mom();
    x.momErr_ = ki.momerr();
    x.dEdgeU_ = dEdgeU;
    x.dEdgeV_ = dEdgeV;
    x.edgeSig_ = (x.sigU() > 0 && x.sigV() > 0) ? std::min(x.dEdgeU_ / x.sigU(), x.dEdgeV_ / x.sigV()) : 0.f;
    if(layerProbe_) probe(x, speed);
    return true;
  }

  TrkCrvMatch TrkCrvMatchTool::match(art::Ptr<KalSeed> const& seed, std::vector<art::Ptr<CrvCoincidenceCluster>> const& usedInFit) const {
    auto const& ks = *seed;
    TrkCrvMatch m;
    m.track_ = seed;
    double t0 = ks.t0Val();
    double t0err = std::sqrt(std::max(ks.t0Var(), 0.0));
    // reference time: the earliest intersection on any reference surface (tracker entry), or the first surface in
    // list order present, nearest in time to t0; KalSeed t0 if none
    int kref = -1;
    if(earliestRef_) {
      for(size_t k = 0; k < ks.intersections().size(); ++k) {
        auto const& ki = ks.intersections()[k];
        if(std::find(refSurfaces_.begin(), refSurfaces_.end(), ki.surfaceId()) == refSurfaces_.end()) continue;
        if(kref < 0 || ki.time() < ks.intersections()[kref].time()) kref = static_cast<int>(k);
      }
    } else {
      for(auto const& rs : refSurfaces_) {
        double bdt = 1e30;
        for(size_t k = 0; k < ks.intersections().size(); ++k) {
          auto const& ki = ks.intersections()[k];
          if(ki.surfaceId() != rs || std::fabs(ki.time() - t0) >= bdt) continue;
          bdt = std::fabs(ki.time() - t0);
          kref = static_cast<int>(k);
        }
        if(kref >= 0) break;
      }
    }
    m.tRef_ = t0;
    bool hasRefPos = kref >= 0;
    KinKal::VEC3 refPos;
    if(hasRefPos) {
      auto const& ki = ks.intersections()[kref];
      m.refSurface_ = ki.surfaceId();
      m.tRef_ = ki.time();
      refPos = ki.position3();
    }
    m.timeOffset_ = timeOffset_;
    if(layerProbe_ && doGeometric_)
      m.probe_ = !pulses_ ? TrkCrvMatch::probeNoPulses : (probeAllPulses_ ? TrkCrvMatch::probeAllPulses : TrkCrvMatch::probeClusteredPulses);

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
        x.cov_ = errorModel(trkCrvPlaneCovariance(ki.pstate_.stateCovariance(), ki.momentum3(), g.u, g.v, g.w, 0.0));
        x.sigT_ = sigT;
        x.mom_ = ki.mom();
        x.momErr_ = ki.momerr();
        KinKal::VEC3 loc = p - g.center;
        x.dEdgeU_ = g.uHalf - std::fabs(loc.Dot(g.u));
        x.dEdgeV_ = g.vHalf - std::fabs(loc.Dot(g.v));
        x.edgeSig_ = (x.sigU() > 0 && x.sigV() > 0) ? std::min(x.dEdgeU_ / x.sigU(), x.dEdgeV_ / x.sigV()) : 0.f;
        if(layerProbe_) probe(x, std::sqrt(ki.velocity().Mag2()));
        m.xings_.push_back(x);
      }
    }
    // a straight-line fit with no CRV crossing passes outside every sector: project it onto the extended plane of each
    // cluster's sector, so that edge clusters get geometric candidates
    bool projectLine = doGeometric_ && m.xings_.empty() && ks.kinematicLineFit();
    std::map<int,int16_t> projected;   // cluster shield -> its projected crossing in m.xings_
    if(projectLine) {
      std::set<int> shields;
      for(auto const& ci : cinfo_)
        if(ci.shield >= 0 && static_cast<size_t>(ci.shield) < sectors_.size() && sectors_[ci.shield].valid) shields.insert(ci.shield);
      for(int s : shields) {
        TrkCrvXing x;
        if(!projectXing(ks, s, sigT, x)) continue;
        projected[s] = static_cast<int16_t>(m.xings_.size());
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
      c.dtRefToF_ = hasRefPos ? c.dtRef_ - std::sqrt((refPos - ci.pos).Mag2()) / CLHEP::c_light : c.dtRef_;
      if(ci.singleEnded) c.setFlag(TrkCrvCandidate::singleEnded);
      if(ci.multiSector) c.setFlag(TrkCrvCandidate::multiSector);
      bool timeOK = keepAll_ || (c.dtRef_ >= dtRefLo_ && c.dtRef_ <= dtRefHi_);
      if(timeOK) c.setFlag(TrkCrvCandidate::timeCand);
      // geometric evidence: project from the CRV crossing nearest in path to the cluster's plane (straight line outside the field)
      int bestx = -1; double bds = 1e30;
      bool hasPlane = ci.shield >= 0 && static_cast<size_t>(ci.shield) < sectors_.size() && sectors_[ci.shield].valid;
      if(hasPlane && projectLine) {   // a projected line: from the projection onto the cluster's own sector
        auto ip = projected.find(ci.shield);
        if(ip != projected.end()) bestx = ip->second;
      } else if(hasPlane) {
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
        // track covariance on the cluster's plane: the reference crossing's covariance carried over the straight step;
        // a projected crossing carries that of the intersection it was projected from, over the step from there
        double dsCov = ki.surfaceId() == x.sid_ ? ds : (ci.pos - ki.position3()).Dot(g.w) / dn;
        c.cov_ = errorModel(trkCrvPlaneCovariance(ki.pstate_.stateCovariance(), ki.momentum3(), g.u, g.v, g.w, dsCov));
        double sv = resolution(sigmaV_, ci.sectorType), su = resolution(sigmaU_, ci.sectorType);
        double st = ci.singleEnded ? sigmaTSingle_ : resolution(sigmaT_, ci.sectorType);
        c.chi2_ = c.dtXing_*c.dtXing_/(st*st + x.sigT_*x.sigT_);
        if(ci.singleEnded) {   // no along-bar measurement: across-bar term only
          c.chi2_ += c.dv_*c.dv_/(sv*sv + c.cov_(1,1));
          c.ndof_ = 2;
        } else {               // full (u, v) covariance plus the CRV resolution on the diagonal
          double a = c.cov_(0,0) + su*su, b = c.cov_(0,1), d = c.cov_(1,1) + sv*sv, det = a*d - b*b;
          c.chi2_ += det > 0 ? (d*c.du_*c.du_ - 2*b*c.du_*c.dv_ + a*c.dv_*c.dv_)/det : 0.0;
          c.ndof_ = 3;
        }
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
      if(!c.hasFlag(TrkCrvCandidate::geomCand)) { c.chi2_ = -1; c.ndof_ = 0; c.du_ = c.dv_ = c.dtXing_ = 0; c.cov_ = TrkCrvPlaneCov(); }
      if(c.hasFlag(TrkCrvCandidate::geomCand) || timeOK) {
        for(auto const& u : usedInFit) if(sameCluster(u, c.cluster_)) { c.setFlag(TrkCrvCandidate::usedInFit); break; }
        m.cands_.push_back(c);
      }
    }
    if(!projected.empty()) {   // keep only the projected crossings a geometric candidate refers to, in sector order
      std::vector<bool> used(m.xings_.size(), false);
      for(auto const& c : m.cands_) if(c.xing_ >= 0) used[c.xing_] = true;
      std::vector<int16_t> remap(m.xings_.size(), -1);
      std::vector<TrkCrvXing> kept;
      for(size_t k = 0; k < m.xings_.size(); ++k)
        if(used[k]) { remap[k] = static_cast<int16_t>(kept.size()); kept.push_back(std::move(m.xings_[k])); }
      m.xings_ = std::move(kept);
      for(auto& c : m.cands_) if(c.xing_ >= 0) c.xing_ = remap[c.xing_];
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
    std::map<std::pair<size_t,int>,bool> usedTrackSector;   // (track, cluster sector): one cluster per track per sector
    for(auto const& [chi2, r] : geomPairs) {
      auto& c = out[r.trk].cands_[r.cand];
      auto ts = std::make_pair(r.trk, static_cast<int>(c.sector_));
      if(usedCluster[c.cluster_.key()] || usedTrackSector[ts]) continue;
      usedCluster[c.cluster_.key()] = usedTrackSector[ts] = true;
      c.setFlag(TrkCrvCandidate::exclusive);
    }
  }
}
