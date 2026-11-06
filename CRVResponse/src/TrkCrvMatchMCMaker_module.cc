//
//  MC truth for TrkCrvMatch by SimParticle identity. For each track: its particle (KalSeedMC primary contributor), the
//  true CRV sector crossings of that particle and its ancestors (MCTrajectory chords through the sector mid-plane
//  rectangles, gap-inclusive, plus CrvStep deposits), the truth of each reco crossing, and the relation of each
//  candidate cluster's deposits to the track's particle. The geometry service must hold the SIMULATION geometry.
//  Particles are identified by (simStage, key), so steps and trajectories from different compressed SimParticle
//  collections of one stage compare correctly.
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/CRVReco/inc/TrkCrvPlaneError.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/MCDataProducts/inc/CrvCoincidenceClusterMCAssns.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/KalSeedMC.hh"
#include "Offline/MCDataProducts/inc/MCTrajectory.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/TrkCrvMatchMC.hh"
#include "Offline/RecoDataProducts/inc/TrkCrvMatch.hh"
#include "KinKal/General/Vectors.hh"

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
      struct SectorGeom { bool valid = false; std::shared_ptr<KinKal::Rectangle> rect; SurfaceId sid; };
      struct Counts {           // per shield
        unsigned recoX = 0, recoTrueDep = 0, recoTrueGap = 0, recoTrueAnc = 0, recoFake = 0, recoFakeEdge = 0;
        unsigned recoFakeNbr = 0, trueMissedNbr = 0;
        unsigned trueClu = 0, trueCluGeom = 0, trueCluExcl = 0;   // the particle's own cluster in this sector
        unsigned trueCluBest = 0, trueCluProj = 0;                // ... best in its sector; projected from another sector
        double chi2nd = 0, chi2ndProj = 0;                        // sums of chi2/ndof of those geometric candidates
        unsigned trueX = 0, trueFound = 0, trueMissed = 0, trueNoGeom = 0;
      };
      // relation of particle q to the track's particle: walk q's ancestry
      uint8_t relation(art::Ptr<SimParticle> q, std::set<PKey> const& self, std::set<PKey> const& anc, PKey const& trkRoot) const;
      PKey root(art::Ptr<SimParticle> p) const;

      art::InputTag matchTag_, ksmcTag_, spTag_, trajTag_, stepTag_, ccmcTag_;
      double edgeSigCut_;
      unsigned maxGen_;
      bool printSummary_;
      int diag_;
      std::vector<SectorGeom> sectors_;  // indexed by shield
      std::vector<int> planeGroup_;      // per shield: sectors sharing one plane get one id, -1 = no rectangle
      std::vector<std::string> shieldNames_;
      std::map<int,Counts> counts_;
      std::array<unsigned,6> bestGeomRel_{}, bestTimeRel_{};
      unsigned nTracks_ = 0, nTracksMC_ = 0, nTracksCrossing_ = 0;
  };

  TrkCrvMatchMCMaker::TrkCrvMatchMCMaker(const Parameters& conf) :
    art::EDProducer(conf),
    matchTag_(conf().matchTag()), ksmcTag_(conf().kalSeedMCAssnsTag()), spTag_(conf().simParticleTag()),
    trajTag_(conf().mcTrajectoryTag()), stepTag_(conf().crvStepTag()), ccmcTag_(conf().clusterMCAssnsTag()),
    edgeSigCut_(conf().edgeSigCut()), maxGen_(conf().maxGenerations()), printSummary_(conf().printSummary()), diag_(conf().diag())
  {
    consumes<TrkCrvMatchCollection>(matchTag_);
    consumes<KalSeedMCAssns>(ksmcTag_);
    consumes<SimParticleCollection>(spTag_);
    consumes<MCTrajectoryCollection>(trajTag_);
    consumes<CrvStepCollection>(stepTag_);
    consumes<CrvCoincidenceClusterMCAssns>(ccmcTag_);
    produces<TrkCrvMatchMCCollection>();
  }

  void TrkCrvMatchMCMaker::beginRun(art::Run&) {
    GeomHandle<CosmicRayShield> crs;
    GeomHandle<KinKalGeom> kkg;
    auto const& shields = crs->getCRSScintillatorShields();
    sectors_.assign(shields.size(), SectorGeom());
    shieldNames_.clear();
    std::map<std::string,int> nameToShield;
    for(size_t i = 0; i < shields.size(); ++i) { nameToShield[shields[i].getName()] = static_cast<int>(i); shieldNames_.push_back(shields[i].getName()); }
    for(auto const& sec : kkg->CRV()->sectors()) {
      auto it = nameToShield.find(sec.sname_);
      if(it == nameToShield.end()) continue;
      auto& g = sectors_[it->second];
      g.valid = true;
      g.rect = sec.sector_;
      try { g.sid = SurfaceId(sec.sname_); } catch(std::exception const&) {}
    }
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

  TrkCrvMatchMCMaker::PKey TrkCrvMatchMCMaker::root(art::Ptr<SimParticle> p) const {
    for(int i = 0; i < 1000 && p->parent().isNonnull(); ++i) p = p->parent();
    return pkey(p);
  }

  uint8_t TrkCrvMatchMCMaker::relation(art::Ptr<SimParticle> q, std::set<PKey> const& self, std::set<PKey> const& anc, PKey const& trkRoot) const {
    if(q.isNull()) return TrkCrvTruth::Background;
    bool family = false;
    for(int i = 0; i < 1000 && q.isNonnull(); ++i) {
      auto k = pkey(q);
      if(self.count(k)) return TrkCrvTruth::SameParticle; // the track's particle (any stage) or a descendant of it
      if(anc.count(k)) family = true;                    // an ancestor or another descendant of an ancestor
      if(q->parent().isNull()) return family ? TrkCrvTruth::SameFamily : (k == trkRoot ? TrkCrvTruth::SamePrimary : TrkCrvTruth::Unrelated);
      q = q->parent();
    }
    return family ? TrkCrvTruth::SameFamily : TrkCrvTruth::Unrelated;
  }

  void TrkCrvMatchMCMaker::produce(art::Event& event) {
    auto const& mH = event.getValidHandle<TrkCrvMatchCollection>(matchTag_);
    auto const& ksmcH = event.getValidHandle<KalSeedMCAssns>(ksmcTag_);
    auto spH = event.getHandle<SimParticleCollection>(spTag_);
    auto trajH = event.getHandle<MCTrajectoryCollection>(trajTag_);
    auto stepH = event.getHandle<CrvStepCollection>(stepTag_);
    auto ccmcH = event.getHandle<CrvCoincidenceClusterMCAssns>(ccmcTag_);
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
    if(ccmcH.isValid()) for(auto const& a : *ccmcH) {
      ccmc[a.first.key()] = a.second;
      clusterPid = a.first.id();
      std::map<int,double> pes;
      for(auto const& rp : a.first->GetCrvRecoPulses()) pes[crs->getBar(rp->GetScintillatorBarIndex()).id().getShieldNumber()] += rp->GetPEs();
      double best = -1;
      for(auto const& [s, pe] : pes) if(pe > best) { best = pe; clusterShield[a.first.key()] = s; }
    }

    auto out = std::make_unique<TrkCrvMatchMCCollection>();
    for(size_t imatch = 0; imatch < mH->size(); ++imatch) {
      auto const& m = (*mH)[imatch];
      TrkCrvMatchMC mc;
      mc.match_ = art::Ptr<TrkCrvMatch>(mH, imatch);
      mc.recoXingTruth_.assign(m.xings_.size(), -1);
      mc.recoXingNeighbor_.assign(m.xings_.size(), -1);
      mc.cands_.assign(m.cands_.size(), TrkCrvCandidateMC());
      ++nTracks_;
      auto im = ksmc.find({m.track_.id(), m.track_.key()});
      bool ok = im != ksmc.end() && spH.isValid();
      if(ok && !im->second->simParticles().empty()) mc.trackSim_ = im->second->simParticle().simParticle(spH);
      if(mc.trackSim_.isNull() || !mc.trackSim_.isAvailable()) { out->push_back(std::move(mc)); continue; }
      ++nTracksMC_;
      PKey trk = pkey(mc.trackSim_);
      PKey trkRoot = root(mc.trackSim_);

      // lineage: generation 0 = the track's particle (and its earlier-stage selves), then ancestors
      std::vector<std::pair<art::Ptr<SimParticle>,uint8_t>> lineage;
      std::set<PKey> anc, self;
      art::Ptr<SimParticle> p = mc.trackSim_;
      uint8_t gen = 0;
      for(int i = 0; i < 1000 && p.isNonnull() && gen <= maxGen_; ++i) {
        lineage.push_back({p, gen});
        if(gen > 0) anc.insert(pkey(p)); else self.insert(pkey(p));
        bool self = p->selfParent();
        p = p->parent();
        if(!self) ++gen;
      }

      // true crossings: trajectory chords through the sector rectangles, then deposits
      std::map<std::pair<size_t,int>, TrkCrvXingMC> xmap;    // (lineage index, shield)
      for(size_t il = 0; il < lineage.size(); ++il) {
        auto const& [sp, g] = lineage[il];
        if(trajH.isValid()) {
          auto it = trajH->find(sp);
          if(it != trajH->end()) {
            auto const& pts = it->second.points();
            double mass = sp->startMomentum().m();
            for(size_t j = 1; j < pts.size(); ++j) {
              KinKal::VEC3 a = toDet(pts[j-1].pos()), b = toDet(pts[j].pos());
              for(size_t s = 0; s < sectors_.size(); ++s) {
                auto const& sg = sectors_[s];
                if(!sg.valid) continue;
                double da = (a - sg.rect->center()).Dot(sg.rect->normal()), db = (b - sg.rect->center()).Dot(sg.rect->normal());
                if(da * db > 0 || da == db) continue;
                double f = da / (da - db);
                KinKal::VEC3 x = a + f * (b - a);
                if(!sg.rect->inBounds(x, 0.0)) continue;
                auto key = std::make_pair(il, static_cast<int>(s));
                double t = pts[j-1].t() + f * (pts[j].t() - pts[j-1].t());
                auto& xm = xmap[key];
                if(xm.geom_ && xm.time_ <= t) continue;    // keep the earliest crossing of a sector
                xm.geom_ = true;
                xm.pos_ = XYZVectorF(x.X(), x.Y(), x.Z());
                auto d = (b - a).Unit();
                xm.dir_ = XYZVectorF(d.X(), d.Y(), d.Z());
                xm.time_ = t;
                double ke = pts[j-1].kineticEnergy() + f * (pts[j].kineticEnergy() - pts[j-1].kineticEnergy());
                xm.mom_ = std::sqrt(std::max(ke * (ke + 2.0 * mass), 0.0));
                auto loc = x - sg.rect->center();
                xm.dEdgeU_ = sg.rect->uHalfLength() - std::fabs(loc.Dot(sg.rect->uDirection()));
                xm.dEdgeV_ = sg.rect->vHalfLength() - std::fabs(loc.Dot(sg.rect->vDirection()));
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
        xm.generation_ = lineage[key.first].second;
        xm.sim_ = lineage[key.first].first;
        if(key.second >= 0 && static_cast<size_t>(key.second) < sectors_.size()) xm.sid_ = sectors_[key.second].sid;
        mc.xings_.push_back(xm);
      }

      // the cluster with the most energy from each crossing's particle in that sector
      if(ccmcH.isValid()) {
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
            if(static_cast<int>(m.cands_[ic].cluster_.key()) == xm.trueCluster_) { xm.trueCand_ = static_cast<int16_t>(ic); break; }
        }
      }

      // reco crossings <-> true crossings by sector; prefer the track's own particle, then the nearest ancestor
      for(size_t ir = 0; ir < m.xings_.size(); ++ir) {
        int best = -1;
        for(size_t ix = 0; ix < mc.xings_.size(); ++ix) {
          auto const& xm = mc.xings_[ix];
          if(xm.sector_ != m.xings_[ir].sector_ || !(xm.geom_ || xm.edep_ > 0)) continue;
          if(best < 0 || xm.generation_ < mc.xings_[best].generation_) best = static_cast<int>(ix);
        }
        mc.recoXingTruth_[ir] = static_cast<int16_t>(best);
        if(best >= 0 && mc.xings_[best].recoXing_ < 0) mc.xings_[best].recoXing_ = static_cast<int16_t>(ir);
        if(best < 0) {   // the track's particle crossed a coplanar neighbour instead
          int rs = m.xings_[ir].sector_;
          for(size_t ix = 0; ix < mc.xings_.size(); ++ix) {
            auto const& xm = mc.xings_[ix];
            if(xm.generation_ != 0 || !(xm.geom_ || xm.edep_ > 0) || xm.sector_ == rs) continue;
            if(rs >= 0 && xm.sector_ >= 0 && planeGroup_[rs] >= 0 && planeGroup_[rs] == planeGroup_[xm.sector_]) { mc.recoXingNeighbor_[ir] = static_cast<int16_t>(ix); break; }
          }
        }
      }

      // candidate clusters: relation of their deposits to the track's particle
      if(ccmcH.isValid() && !m.cands_.empty() && m.cands_.front().cluster_.id() != clusterPid && !ccmc.empty())
        mf::LogWarning("TrkCrvMatchMCMaker") << "cluster MC associations are over a different cluster collection than the matcher used";
      for(size_t ic = 0; ic < m.cands_.size(); ++ic) {
        auto& cm = mc.cands_[ic];
        auto it = ccmc.find(m.cands_[ic].cluster_.key());
        if(it == ccmc.end() || m.cands_[ic].cluster_.id() != clusterPid) continue;
        if(!it->second->HasMCInfo()) { cm.rel_ = TrkCrvTruth::Background; continue; }
        uint8_t rel = TrkCrvTruth::Background;
        for(auto const& pi : it->second->GetPulses()) {
          if(pi._energyDeposited <= 0) continue;
          uint8_t r = relation(pi._simParticle, self, anc, trkRoot);
          cm.edepTotal_ += pi._energyDeposited;
          if(r == TrkCrvTruth::SameParticle) cm.edepSame_ += pi._energyDeposited;
          if(r == TrkCrvTruth::SameFamily) cm.edepFamily_ += pi._energyDeposited;
          rel = std::min(rel, r);
        }
        cm.rel_ = rel;
      }

      // bookkeeping for the end-of-job summary
      bool crosses = false;
      for(size_t ir = 0; ir < m.xings_.size(); ++ir) {
        auto& c = counts_[m.xings_[ir].sector_];
        ++c.recoX;
        int t = mc.recoXingTruth_[ir];
        if(t < 0) { ++c.recoFake; if(m.xings_[ir].edgeSig_ < edgeSigCut_) ++c.recoFakeEdge; if(mc.recoXingNeighbor_[ir] >= 0) ++c.recoFakeNbr; }
        else if(mc.xings_[t].generation_ > 0) ++c.recoTrueAnc;
        else if(mc.xings_[t].edep_ > 0) ++c.recoTrueDep;
        else ++c.recoTrueGap;
      }
      for(auto const& xm : mc.xings_) {
        if(xm.generation_ != 0) continue;
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
        if(xm.generation_ != 0 || xm.trueCluster_ < 0) continue;
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
      if(m.bestGeom_ >= 0) ++bestGeomRel_[mc.cands_[m.bestGeom_].rel_];
      if(m.bestTime_ >= 0) ++bestTimeRel_[mc.cands_[m.bestTime_].rel_];

      if(diag_ > 0) {
        mf::LogInfo log("TrkCrvMatchMCMaker");
        log << "event " << event.id() << " track " << m.track_.key() << " sim " << mc.trackSim_->pdgId() << " stage/key "
            << trk.first << "/" << trk.second << " reco xings " << m.xings_.size() << " true xings " << mc.xings_.size();
        for(size_t ir = 0; ir < m.xings_.size(); ++ir)
          log << "\n  reco " << m.xings_[ir].sid_ << " inb " << m.xings_[ir].inBounds_ << " edgeU " << m.xings_[ir].dEdgeU_
              << " edgeV " << m.xings_[ir].dEdgeV_ << " edgeSig " << m.xings_[ir].edgeSig_ << " truth " << mc.recoXingTruth_[ir] << " nbr " << mc.recoXingNeighbor_[ir];
        for(auto const& xm : mc.xings_)
          log << "\n  true shield " << xm.sector_ << " gen " << int(xm.generation_) << " geom " << xm.geom_ << " edep " << xm.edep_
              << " layers " << int(xm.nLayers_) << " edgeU " << xm.dEdgeU_ << " edgeV " << xm.dEdgeV_
              << " recoXing " << xm.recoXing_ << " cluster " << xm.trueCluster_;
        // pulls of own-cluster geometric candidates: matcher residuals, and the raw-error track against the true crossing
        for(auto const& xm : mc.xings_) {
          if(xm.generation_ != 0 || !xm.geom_ || xm.trueCand_ < 0) continue;
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
          log << "\n  pull shield " << xm.sector_ << " same " << c.hasFlag(TrkCrvCandidate::sameSector)
              << " single " << c.hasFlag(TrkCrvCandidate::singleEnded) << " s " << s << " dw " << dw
              << " du " << c.du_ << " dv " << c.dv_ << " sigU " << c.sigU_ << " sigV " << c.sigV_
              << " tdu " << r.Dot(rect.uDirection()) << " tdv " << r.Dot(rect.vDirection())
              << " teu " << trkCrvPlaneError(cov, mom, rect.uDirection(), w, s)
              << " tev " << trkCrvPlaneError(cov, mom, rect.vDirection(), w, s)
              << " teu0 " << trkCrvPlaneError(cov, mom, rect.uDirection(), w, 0.0)
              << " tev0 " << trkCrvPlaneError(cov, mom, rect.vDirection(), w, 0.0);
        }
      }
      out->push_back(std::move(mc));
    }
    event.put(std::move(out));
  }

  void TrkCrvMatchMCMaker::endJob() {
    if(!printSummary_) return;
    static const char* relName[6] = {"SameParticle","SameFamily","SamePrimary","Unrelated","Background","NoMC"};
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
    for(int i = 0; i < 6; ++i) log << " " << relName[i] << " " << bestGeomRel_[i];
    log << "\nbest time candidate:     ";
    for(int i = 0; i < 6; ++i) log << " " << relName[i] << " " << bestTimeRel_[i];
  }
}

DEFINE_ART_MODULE(mu2e::TrkCrvMatchMCMaker)
