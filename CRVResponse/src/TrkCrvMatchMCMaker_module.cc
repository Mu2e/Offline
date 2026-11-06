//
//  MC truth for TrkCrvMatch. For each track: its KalSeedMC and particle (primary contributor), the true CRV sector
//  crossings of that particle and its ancestors (MCTrajectory segments through the sector mid-plane rectangles,
//  gap-inclusive, plus CrvStep deposits), the truth of each reco crossing, and the MCRelationship of each candidate
//  cluster's deposits to the track's particle. The geometry service must hold the SIMULATION geometry.
//  Steps and trajectories are matched to particles by (simStage, key), so different compressed SimParticle
//  collections of one stage compare correctly.
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/CRVReco/inc/TrkCrvGeom.hh"
#include "Offline/CRVReco/inc/TrkCrvPlaneError.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/MCDataProducts/inc/CrvCoincidenceClusterMCAssns.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/KalSeedMC.hh"
#include "Offline/MCDataProducts/inc/MCTrajectory.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/TrkCrvMatchMC.hh"
#include "Offline/RecoDataProducts/inc/TrkCrvMatch.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Geometry/LineSegment.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace mu2e {

  // relation of an ancestor to the track's particle, from its generation in the stage-aware lineage (MCRelationship
  // from the Ptrs would count each simulation-stage copy as a generation)
  inline MCRelationship ancestorRelation(unsigned gen) {
    if(gen == 0) return MCRelationship(MCRelationship::same, 0);
    return MCRelationship(gen == 1 ? MCRelationship::mother : MCRelationship::umother, static_cast<int8_t>(gen));
  }

  class TrkCrvMatchMCMaker : public art::EDProducer {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> matchTag{Name("TrkCrvMatchCollection"), Comment("TrkCrvMatch collection to add truth to")};
        fhicl::Atom<art::InputTag> kalSeedMCAssnsTag{Name("KalSeedMCAssns"), Comment("KalSeed to KalSeedMC associations"), art::InputTag("SelectReco")};
        fhicl::Atom<art::InputTag> simParticleTag{Name("SimParticleCollection"), Comment("SimParticles the KalSeedMC keys refer to"), art::InputTag("compressRecoMCs")};
        fhicl::Atom<art::InputTag> mcTrajectoryTag{Name("MCTrajectoryCollection"), Comment("MC trajectories for the gap-inclusive crossings"), art::InputTag("compressRecoMCs")};
        fhicl::Atom<art::InputTag> crvStepTag{Name("CrvStepCollection"), Comment("CRV steps for the deposits; the digi-stage collection holds every digitised step"), art::InputTag("compressDigiMCs")};
        fhicl::Atom<bool> requireCrvSteps{Name("RequireCrvSteps"), Comment("throw if the CRV steps are missing; false fills no deposits (hasDeposits_ false), e.g. on reco output, which drops them"), true};
        fhicl::Atom<art::InputTag> clusterMCAssnsTag{Name("CrvCoincidenceClusterMCAssns"), Comment("cluster to cluster-MC associations, over the clusters the matcher used"), art::InputTag("CrvCoincidenceClusterMCAssns")};
        fhicl::Atom<double> edgeSigCut{Name("EdgeSigCut"), Comment("summary: a fake counts as near-edge when its edgeSig is below this"), 3.0};
        fhicl::Atom<unsigned> maxGenerations{Name("MaxGenerations"), Comment("ancestors of the track's particle searched for CRV crossings"), 3};
        fhicl::Atom<bool> printSummary{Name("PrintSummary"), Comment("print per-sector crossing truth counts at end of job"), true};
        fhicl::Atom<int> diag{Name("diagLevel"), Comment("diagnostic printout level"), 0};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit TrkCrvMatchMCMaker(const Parameters& conf);
      void beginRun(art::Run& run) override;
      void produce(art::Event& event) override;
      void endJob() override;

    private:
      using PKey = std::pair<unsigned, unsigned long>;   // (simStage, key): particle identity across collections
      static PKey pkey(art::Ptr<SimParticle> const& p) { return {p->simStage(), p.key()}; }
      static constexpr int nRel = 9;                     // summary bins: no MC, then MCRelationship none ... usibling
      struct Counts {           // per shield
        unsigned recoX = 0, recoTrueDep = 0, recoTrueGap = 0, recoTrueAnc = 0, recoFake = 0, recoFakeEdge = 0;
        unsigned recoFakeNbr = 0, trueMissedNbr = 0;
        unsigned trueClu = 0, trueCluGeom = 0, trueCluExcl = 0;   // the particle's own cluster in this sector
        unsigned trueCluBest = 0, trueCluProj = 0;                // ... best in its sector; projected from another sector
        double chi2nd = 0, chi2ndProj = 0;                        // sums of chi2/ndof of those geometric candidates
        unsigned trueX = 0, trueFound = 0, trueMissed = 0, trueNoGeom = 0;
      };
      art::InputTag matchTag_, ksmcTag_, spTag_, trajTag_, stepTag_, ccmcTag_;
      bool requireSteps_;
      double edgeSigCut_;
      unsigned maxGen_;
      bool printSummary_;
      int diag_;
      std::vector<TrkCrvSector> sectors_;  // indexed by shield
      std::vector<int> planeGroup_;      // per shield: sectors sharing one plane get one id, -1 = no rectangle
      std::vector<std::string> shieldNames_;
      std::map<int,Counts> counts_;
      std::array<unsigned,nRel> bestGeomRel_{}, bestTimeRel_{};
      std::vector<float> dtRefSame_;     // dtRef of time candidates whose cluster is the track's own particle (DtRefCenter)
      unsigned nTracks_ = 0, nTracksMC_ = 0, nTracksCrossing_ = 0;
  };

  TrkCrvMatchMCMaker::TrkCrvMatchMCMaker(const Parameters& conf) :
    art::EDProducer(conf),
    matchTag_(conf().matchTag()), ksmcTag_(conf().kalSeedMCAssnsTag()), spTag_(conf().simParticleTag()),
    trajTag_(conf().mcTrajectoryTag()), stepTag_(conf().crvStepTag()), ccmcTag_(conf().clusterMCAssnsTag()),
    requireSteps_(conf().requireCrvSteps()), edgeSigCut_(conf().edgeSigCut()), maxGen_(conf().maxGenerations()), printSummary_(conf().printSummary()), diag_(conf().diag())
  {
    consumes<TrkCrvMatchCollection>(matchTag_);
    consumes<KalSeedMCAssns>(ksmcTag_);
    consumes<SimParticleCollection>(spTag_);
    consumes<MCTrajectoryCollection>(trajTag_);
    if(requireSteps_) consumes<CrvStepCollection>(stepTag_);
    else mayConsume<CrvStepCollection>(stepTag_);
    consumes<CrvCoincidenceClusterMCAssns>(ccmcTag_);
    produces<TrkCrvMatchMCCollection>();
  }

  void TrkCrvMatchMCMaker::beginRun(art::Run&) {
    sectors_ = trkCrvSectors();
    shieldNames_.clear();
    for(auto const& sg : sectors_) shieldNames_.push_back(sg.name);
    // coplanar sectors: parallel normals and centres within 5 mm of one plane
    planeGroup_.assign(sectors_.size(), -1);
    int ng = 0;
    for(size_t s = 0; s < sectors_.size(); ++s) {
      if(!sectors_[s].valid || planeGroup_[s] >= 0) continue;
      planeGroup_[s] = ng;
      auto const& rs = *sectors_[s].rect;
      for(size_t t = s + 1; t < sectors_.size(); ++t) {
        if(!sectors_[t].valid || planeGroup_[t] >= 0) continue;
        auto const& rt = *sectors_[t].rect;
        if(std::fabs(rs.normal().Dot(rt.normal())) > 0.999 && std::fabs((rt.center() - rs.center()).Dot(rs.normal())) < 5.0) planeGroup_[t] = ng;
      }
      ++ng;
    }
  }

  void TrkCrvMatchMCMaker::produce(art::Event& event) {
    auto const& mH = event.getValidHandle<TrkCrvMatchCollection>(matchTag_);
    auto const& ksmcH = event.getValidHandle<KalSeedMCAssns>(ksmcTag_);
    auto spH = event.getHandle<SimParticleCollection>(spTag_);   // a plain Handle, as SimPartStub::simParticle takes
    if(!spH.isValid()) throw cet::exception("TrkCrvMatchMCMaker") << "no SimParticleCollection " << spTag_;
    auto const& trajH = event.getValidHandle<MCTrajectoryCollection>(trajTag_);
    auto const& ccmcH = event.getValidHandle<CrvCoincidenceClusterMCAssns>(ccmcTag_);
    auto stepH = event.getHandle<CrvStepCollection>(stepTag_);
    if(!stepH.isValid() && requireSteps_)
      throw cet::exception("TrkCrvMatchMCMaker") << "no CrvStepCollection " << stepTag_ << " (reco output drops it; set RequireCrvSteps : false to run without deposits)";
    GeomHandle<CosmicRayShield> crs;
    GeomHandle<DetectorSystem> det;
    auto toDet = [&](CLHEP::Hep3Vector const& v) { auto d = det->toDetector(v); return KinKal::VEC3(d.x(), d.y(), d.z()); };

    // (KalSeed collection, key) -> KalSeedMC: SelectRecoMC writes one Assns over every KalSeed collection it selects
    std::map<std::pair<art::ProductID, size_t>, art::Ptr<KalSeedMC>> ksmc;
    for(auto const& a : *ksmcH) ksmc[{a.first.id(), a.first.key()}] = a.second;
    // cluster key -> cluster MC, and each cluster's dominant shield
    std::map<size_t, art::Ptr<CrvCoincidenceClusterMC>> ccmc;
    std::map<size_t, int> clusterShield;
    art::ProductID clusterPid;
    for(auto const& a : *ccmcH) {
      ccmc[a.first.key()] = a.second;
      clusterPid = a.first.id();
      clusterShield[a.first.key()] = trkCrvDominantShield(*a.first, *crs);
    }
    // a deposit's particle in the track's SimParticle collection (same key), so MCRelationship compares like Ptrs
    auto inTrackCollection = [&](art::Ptr<SimParticle> const& q) {
      return (q.isNonnull() && spH->has(cet::map_vector_key(q.key()))) ? art::Ptr<SimParticle>(spH, q.key()) : art::Ptr<SimParticle>();
    };

    auto out = std::make_unique<TrkCrvMatchMCCollection>();
    for(size_t imatch = 0; imatch < mH->size(); ++imatch) {
      auto const& m = (*mH)[imatch];
      TrkCrvMatchMC mc;
      mc.match_ = art::Ptr<TrkCrvMatch>(mH, imatch);
      mc.recoXingTruth_.assign(m.xings_.size(), -1);
      mc.recoXingNeighbor_.assign(m.xings_.size(), -1);
      mc.cands_.assign(m.cands_.size(), TrkCrvCandidateMC());
      ++nTracks_;
      mc.hasDeposits_ = stepH.isValid();
      auto im = ksmc.find({m.track_.id(), m.track_.key()});
      art::Ptr<SimParticle> trackSim;
      if(im != ksmc.end()) {
        mc.kalSeedMC_ = im->second;
        if(!im->second->simParticles().empty()) trackSim = im->second->simParticle().simParticle(spH);
      }
      if(trackSim.isNull() || !trackSim.isAvailable()) { out->push_back(std::move(mc)); continue; }
      ++nTracksMC_;
      PKey trk = pkey(trackSim);

      // lineage: generation 0 = the track's particle (and its earlier-stage selves), then ancestors
      std::vector<std::pair<art::Ptr<SimParticle>,uint8_t>> lineage;
      art::Ptr<SimParticle> p = trackSim;
      uint8_t gen = 0;
      for(int i = 0; i < 1000 && p.isNonnull() && gen <= maxGen_; ++i) {
        lineage.push_back({p, gen});
        bool selfParent = p->selfParent();
        p = p->parent();
        if(!selfParent) ++gen;
      }

      // true crossings: trajectory segments through the sector rectangles, then deposits
      std::map<std::pair<size_t,int>, TrkCrvXingMC> xmap;    // (lineage index, shield)
      for(size_t il = 0; il < lineage.size(); ++il) {
        auto const& [sp, g] = lineage[il];
        {
          auto it = trajH->find(sp);
          if(it != trajH->end()) {
            auto const& pts = it->second.points();
            double mass = sp->startMomentum().m();
            for(size_t j = 1; j < pts.size(); ++j) {
              KinKal::VEC3 a = toDet(pts[j-1].pos()), b = toDet(pts[j].pos());
              if((b - a).Mag2() <= 0) continue;
              KinKal::LineSegment seg(a, b);   // the trajectory between two saved points
              for(size_t s = 0; s < sectors_.size(); ++s) {
                auto const& sg = sectors_[s];
                if(!sg.valid) continue;
                double dist = 0;
                auto flag = sg.rect->intersect(seg, dist, true, 0.0);
                if(!flag.onsurface_ || !flag.inbounds_ || dist < 0 || dist > seg.length()) continue;
                double f = dist / seg.length();
                KinKal::VEC3 x = seg.position(dist);
                auto key = std::make_pair(il, static_cast<int>(s));
                double t = pts[j-1].t() + f * (pts[j].t() - pts[j-1].t());
                auto& xm = xmap[key];
                if(xm.geom_ && xm.time_ <= t) continue;    // keep the earliest crossing of a sector
                xm.geom_ = true;
                xm.pos_ = XYZVectorF(x.X(), x.Y(), x.Z());
                auto const& d = seg.direction();
                xm.dir_ = XYZVectorF(d.X(), d.Y(), d.Z());
                xm.time_ = t;
                double ke = pts[j-1].kineticEnergy() + f * (pts[j].kineticEnergy() - pts[j-1].kineticEnergy());
                xm.mom_ = std::sqrt(std::max(ke * (ke + 2.0 * mass), 0.0));
                auto loc = x - sg.center;
                xm.dEdgeU_ = sg.uHalf - std::fabs(loc.Dot(sg.u));
                xm.dEdgeV_ = sg.vHalf - std::fabs(loc.Dot(sg.v));
              }
            }
          }
        }
        if(stepH.isValid()) {
          std::map<int, std::set<int>> layers;
          std::map<int, KinKal::VEC3> wpos;
          PKey k = pkey(sp);
          for(auto const& st : *stepH) {
            if(st.simParticle().isNull() || pkey(st.simParticle()) != k || st.visibleEDep() <= 0) continue;
            auto const& bid = crs->getBar(st.barIndex()).id();
            int s = bid.getShieldNumber();
            auto& xm = xmap[std::make_pair(il, s)];
            xm.edep_ += st.visibleEDep();
            layers[s].insert(bid.getLayerNumber());
            wpos[s] += st.visibleEDep() * 0.5 * (toDet(st.startPosition()) + toDet(st.endPosition()));
          }
          for(auto const& [s, ls] : layers) {
            auto& xm = xmap[std::make_pair(il, s)];
            xm.nLayers_ = static_cast<uint8_t>(ls.size());
            auto w = wpos[s] / xm.edep_;
            xm.depPos_ = XYZVectorF(w.X(), w.Y(), w.Z());
          }
        }
      }
      for(auto& [key, xm] : xmap) {
        xm.sector_ = static_cast<int16_t>(key.second);
        xm.rel_ = ancestorRelation(lineage[key.first].second);
        xm.sim_ = lineage[key.first].first;
        if(key.second >= 0 && static_cast<size_t>(key.second) < sectors_.size()) xm.sid_ = sectors_[key.second].sid;
        mc.xings_.push_back(xm);
      }

      // the cluster with the most energy from each crossing's particle in that sector
      {
        for(size_t ix = 0; ix < mc.xings_.size(); ++ix) {
          auto& xm = mc.xings_[ix];
          PKey k = pkey(xm.sim_);
          double best = 0;
          for(auto const& [ckey, cmc] : ccmc) {
            auto cs = clusterShield.find(ckey);
            if(cs == clusterShield.end() || cs->second != xm.sector_ || !cmc->HasMCInfo()) continue;
            double e = 0;
            for(auto const& pi : cmc->GetPulses()) if(pi._simParticle.isNonnull() && pkey(pi._simParticle) == k) e += pi._energyDeposited;
            if(e > best) { best = e; xm.trueCluster_ = static_cast<int16_t>(ckey); }
          }
          for(size_t ic = 0; ic < m.cands_.size() && xm.trueCluster_ >= 0; ++ic)
            if(m.cands_[ic].cluster_.id() == clusterPid && static_cast<int>(m.cands_[ic].cluster_.key()) == xm.trueCluster_) { xm.trueCand_ = static_cast<int16_t>(ic); break; }
        }
      }

      // reco crossings <-> true crossings by sector; prefer the track's own particle, then the nearest ancestor
      for(size_t ir = 0; ir < m.xings_.size(); ++ir) {
        int best = -1;
        for(size_t ix = 0; ix < mc.xings_.size(); ++ix) {
          auto const& xm = mc.xings_[ix];
          if(xm.sector_ != m.xings_[ir].sector_ || !(xm.geom_ || xm.edep_ > 0)) continue;
          if(best < 0 || xm.rel_.removal() < mc.xings_[best].rel_.removal()) best = static_cast<int>(ix);
        }
        mc.recoXingTruth_[ir] = static_cast<int16_t>(best);
        if(best >= 0 && mc.xings_[best].recoXing_ < 0) mc.xings_[best].recoXing_ = static_cast<int16_t>(ir);
        if(best < 0) {   // the track's particle crossed a coplanar neighbour instead
          int rs = m.xings_[ir].sector_;
          for(size_t ix = 0; ix < mc.xings_.size(); ++ix) {
            auto const& xm = mc.xings_[ix];
            if(xm.rel_ != MCRelationship::same || !(xm.geom_ || xm.edep_ > 0) || xm.sector_ == rs) continue;
            if(rs >= 0 && xm.sector_ >= 0 && planeGroup_[rs] >= 0 && planeGroup_[rs] == planeGroup_[xm.sector_]) { mc.recoXingNeighbor_[ir] = static_cast<int16_t>(ix); break; }
          }
        }
      }

      // candidate clusters: the closest MCRelationship of their deposit particles to the track's particle. A particle
      // simulated in several stages has a copy per stage (self-parent chain), and a deposit can point at any of them,
      // so each deposit is compared with every copy of the track's particle and the closest relation is kept
      std::vector<art::Ptr<SimParticle>> trackSelves;
      for(auto const& [sp, g] : lineage) if(g == 0) trackSelves.push_back(sp);
      if(!m.cands_.empty() && m.cands_.front().cluster_.id() != clusterPid && !ccmc.empty())
        mf::LogWarning("TrkCrvMatchMCMaker") << "cluster MC associations are over a different cluster collection than the matcher used";
      for(size_t ic = 0; ic < m.cands_.size(); ++ic) {
        auto& cm = mc.cands_[ic];
        auto it = ccmc.find(m.cands_[ic].cluster_.key());
        if(it == ccmc.end() || m.cands_[ic].cluster_.id() != clusterPid || !it->second->HasMCInfo()) continue;
        cm.hasMC_ = true;
        for(auto const& pi : it->second->GetPulses()) {
          if(pi._energyDeposited <= 0) continue;
          auto q = inTrackCollection(pi._simParticle);
          MCRelationship r;
          for(auto const& self : trackSelves) { MCRelationship rs(q, self); if(r < rs) r = rs; }
          if(diag_ > 1) {
            mf::LogInfo dl("TrkCrvMatchMCMaker");
            dl << "event " << event.id() << " track " << m.track_.key() << " cand " << ic << " deposit sim "
               << (pi._simParticle.isNonnull() ? std::to_string(pi._simParticle.id().value()) + ":" + std::to_string(pi._simParticle.key()) : std::string("null"))
               << " in track collection " << q.isNonnull();
            if(q.isNonnull()) dl << " stage " << q->simStage() << " selfParent " << q->selfParent();
            dl << " rel " << int(r.relationship()) << "/" << int(r.removal()) << " | track copies";
            for(auto const& self : trackSelves) dl << " " << self.id().value() << ":" << self.key() << "(stage " << self->simStage() << ")";
          }
          cm.edepTotal_ += pi._energyDeposited;
          if(r == MCRelationship::same) cm.edepSame_ += pi._energyDeposited;
          if(r != MCRelationship::none) cm.edepRelated_ += pi._energyDeposited;
          if(cm.rel_ < r) cm.rel_ = r;
        }
        if(cm.rel_ == MCRelationship::same && m.cands_[ic].hasFlag(TrkCrvCandidate::timeCand)) dtRefSame_.push_back(m.cands_[ic].dtRef_);
      }

      // bookkeeping for the end-of-job summary
      bool crosses = false;
      for(size_t ir = 0; ir < m.xings_.size(); ++ir) {
        auto& c = counts_[m.xings_[ir].sector_];
        ++c.recoX;
        int t = mc.recoXingTruth_[ir];
        if(t < 0) { ++c.recoFake; if(m.xings_[ir].edgeSig_ < edgeSigCut_) ++c.recoFakeEdge; if(mc.recoXingNeighbor_[ir] >= 0) ++c.recoFakeNbr; }
        else if(mc.xings_[t].rel_ != MCRelationship::same) ++c.recoTrueAnc;
        else if(mc.xings_[t].edep_ > 0) ++c.recoTrueDep;
        else ++c.recoTrueGap;
      }
      for(auto const& xm : mc.xings_) {
        if(xm.rel_ != MCRelationship::same) continue;
        auto& c = counts_[xm.sector_];
        if(!xm.geom_) { if(xm.edep_ > 0) ++c.trueNoGeom; continue; }
        crosses = true;
        ++c.trueX;
        if(xm.recoXing_ >= 0) ++c.trueFound;
        else {
          ++c.trueMissed;
          for(auto const& rx : m.xings_)
            if(rx.sector_ >= 0 && xm.sector_ >= 0 && planeGroup_[rx.sector_] >= 0 && planeGroup_[rx.sector_] == planeGroup_[xm.sector_]) { ++c.trueMissedNbr; break; }
        }
      }
      for(auto const& xm : mc.xings_) {   // is the particle's own cluster in each sector associated geometrically?
        if(xm.rel_ != MCRelationship::same || xm.trueCluster_ < 0) continue;
        auto& c = counts_[xm.sector_];
        ++c.trueClu;
        bool geom = false, excl = false;
        if(xm.trueCand_ >= 0) {
          auto const& cand = m.cands_[xm.trueCand_];
          if(cand.hasFlag(TrkCrvCandidate::geomCand)) {
            geom = true;
            excl = cand.hasFlag(TrkCrvCandidate::exclusive);
            if(cand.hasFlag(TrkCrvCandidate::bestInSector)) ++c.trueCluBest;
            double cn = cand.ndof_ > 0 ? cand.chi2_ / cand.ndof_ : 0.0;
            c.chi2nd += cn;
            if(!cand.hasFlag(TrkCrvCandidate::sameSector)) { ++c.trueCluProj; c.chi2ndProj += cn; }
          }
        }
        if(geom) ++c.trueCluGeom;
        if(excl) ++c.trueCluExcl;
      }
      if(crosses) ++nTracksCrossing_;
      auto relBin = [](TrkCrvCandidateMC const& cm) { return cm.hasMC_ ? cm.rel_.relationship() + 2 : 0; };
      if(m.bestGeom_ >= 0) ++bestGeomRel_[relBin(mc.cands_[m.bestGeom_])];
      if(m.bestTime_ >= 0) ++bestTimeRel_[relBin(mc.cands_[m.bestTime_])];

      if(diag_ > 0) {
        mf::LogInfo log("TrkCrvMatchMCMaker");
        log << "event " << event.id() << " track " << m.track_.key() << " sim " << trackSim->pdgId() << " stage/key "
            << trk.first << "/" << trk.second << " reco xings " << m.xings_.size() << " true xings " << mc.xings_.size();
        for(size_t ir = 0; ir < m.xings_.size(); ++ir)
          log << "\n  reco " << m.xings_[ir].sid_ << " inb " << m.xings_[ir].inBounds_ << " edgeU " << m.xings_[ir].dEdgeU_
              << " edgeV " << m.xings_[ir].dEdgeV_ << " edgeSig " << m.xings_[ir].edgeSig_ << " truth " << mc.recoXingTruth_[ir] << " nbr " << mc.recoXingNeighbor_[ir];
        for(auto const& xm : mc.xings_)
          log << "\n  true shield " << xm.sector_ << " gen " << int(xm.rel_.removal()) << " geom " << xm.geom_ << " edep " << xm.edep_
              << " layers " << int(xm.nLayers_) << " edgeU " << xm.dEdgeU_ << " edgeV " << xm.dEdgeV_
              << " recoXing " << xm.recoXing_ << " cluster " << xm.trueCluster_;
        // pulls of own-cluster geometric candidates: matcher residuals, and the raw-error track against the true crossing
        for(auto const& xm : mc.xings_) {
          if(xm.rel_ != MCRelationship::same || !xm.geom_ || xm.trueCand_ < 0) continue;
          auto const& c = m.cands_[xm.trueCand_];
          if(!c.hasFlag(TrkCrvCandidate::geomCand) || c.xing_ < 0) continue;
          if(xm.sector_ < 0 || static_cast<size_t>(xm.sector_) >= sectors_.size() || !sectors_[xm.sector_].valid) continue;
          auto const& rect = *sectors_[xm.sector_].rect;
          auto const& rx = m.xings_[c.xing_];
          auto const& ki = m.track_->intersections()[rx.inter_];
          KinKal::VEC3 px(rx.pos_.X(), rx.pos_.Y(), rx.pos_.Z());
          KinKal::VEC3 mom = ki.momentum3(), d = mom.Unit();
          KinKal::VEC3 const& w = rect.normal();
          double dw = d.Dot(w);
          if(std::fabs(dw) < 1e-3) continue;
          double s = (rect.center() - px).Dot(w) / dw;
          KinKal::VEC3 r = KinKal::VEC3(xm.pos_.X(), xm.pos_.Y(), xm.pos_.Z()) - (px + s * d);
          auto const& cov = ki.pstate_.stateCovariance();
          auto ce = trkCrvPlaneCovariance(cov, mom, rect.uDirection(), rect.vDirection(), w, s);
          auto c0 = trkCrvPlaneCovariance(cov, mom, rect.uDirection(), rect.vDirection(), w, 0.0);
          log << "\n  pull shield " << xm.sector_ << " same " << c.hasFlag(TrkCrvCandidate::sameSector)
              << " single " << c.hasFlag(TrkCrvCandidate::singleEnded) << " s " << s << " dw " << dw
              << " du " << c.du_ << " dv " << c.dv_ << " sigU " << c.sigU() << " sigV " << c.sigV()
              << " tdu " << r.Dot(rect.uDirection()) << " tdv " << r.Dot(rect.vDirection())
              << " teu " << std::sqrt(ce(0,0)) << " tev " << std::sqrt(ce(1,1)) << " rhouv " << (ce(0,0) > 0 && ce(1,1) > 0 ? ce(0,1)/std::sqrt(ce(0,0)*ce(1,1)) : 0.0)
              << " teu0 " << std::sqrt(c0(0,0)) << " tev0 " << std::sqrt(c0(1,1));
        }
      }
      out->push_back(std::move(mc));
    }
    event.put(std::move(out));
  }

  void TrkCrvMatchMCMaker::endJob() {
    if(!printSummary_) return;
    static const char* relName[nRel] = {"noMC","none","same","daughter","mother","sibling","udaughter","umother","usibling"};
    mf::LogInfo log("TrkCrvMatchMCMaker");
    log << "TrkCrvMatchMC summary (" << moduleDescription().moduleLabel() << "): tracks " << nTracks_ << ", with MC " << nTracksMC_
        << ", whose particle crosses a CRV sector " << nTracksCrossing_ << "\n"
        << std::setw(10) << "sector" << std::setw(8) << "reco" << std::setw(9) << "true:dep" << std::setw(9) << "true:gap"
        << std::setw(9) << "true:anc" << std::setw(7) << "fake" << std::setw(10) << "fake:edge" << std::setw(9) << "fake:nbr" << " |" << std::setw(7) << "truth" << std::setw(7) << "found"
        << std::setw(8) << "missed" << std::setw(9) << "miss:nbr" << std::setw(12) << "dep-no-geom";
    for(auto const& [s, c] : counts_) {
      std::string n = (s >= 0 && static_cast<size_t>(s) < shieldNames_.size()) ? shieldNames_[s] : std::to_string(s);
      log << "\n" << std::setw(10) << n << std::setw(8) << c.recoX << std::setw(9) << c.recoTrueDep << std::setw(9) << c.recoTrueGap
          << std::setw(9) << c.recoTrueAnc << std::setw(7) << c.recoFake << std::setw(10) << c.recoFakeEdge << std::setw(9) << c.recoFakeNbr << " |" << std::setw(7) << c.trueX << std::setw(7) << c.trueFound
          << std::setw(8) << c.trueMissed << std::setw(9) << c.trueMissedNbr << std::setw(12) << c.trueNoGeom;
    }
    for(auto const& [s, c] : counts_) {
      std::string n = (s >= 0 && static_cast<size_t>(s) < shieldNames_.size()) ? shieldNames_[s] : std::to_string(s);
      log << "\nclusters " << n << " true " << c.trueClu << " geom " << c.trueCluGeom << " exclusive " << c.trueCluExcl
          << " best " << c.trueCluBest << " chi2nd " << c.chi2nd << " proj " << c.trueCluProj << " chi2ndProj " << c.chi2ndProj;
    }
    log << "\nbest geometric candidate:";
    for(int i = 0; i < nRel; ++i) log << " " << relName[i] << " " << bestGeomRel_[i];
    log << "\nbest time candidate:     ";
    for(int i = 0; i < nRel; ++i) log << " " << relName[i] << " " << bestTimeRel_[i];
    if(!dtRefSame_.empty()) {
      std::sort(dtRefSame_.begin(), dtRefSame_.end());
      auto q = [&](double f) { return dtRefSame_[std::min(dtRefSame_.size() - 1, static_cast<size_t>(f * dtRefSame_.size()))]; };
      log << "\nown-particle time candidates: " << dtRefSame_.size() << ", dtRef median " << q(0.5) << " ns, 16-84% " << q(0.16) << " " << q(0.84) << " ns";
    }
  }
}

DEFINE_ART_MODULE(mu2e::TrkCrvMatchMCMaker)
