//
// Assign each calorimeter MC energy deposit to its calo-entrant ancestor:
// the highest SimParticle in the Geant4 parent chain that also deposited
// energy in the same calorimeter disk (the shower originator).
//
// For each CaloHitMC in the input collection this module emits a
// CaloHitEntrant whose entrants vector is aligned with
// CaloHitMC::energyDeposits() and which carries an art::Ptr back to its
// CaloHitMC; the output collection is index-parallel to the input
// CaloHitMCCollection.
//
// Disk resolution prefers CaloHitMC::crystalID() when it is filled. For
// files produced before that member existed (reads back -1) the disk is
// resolved through the optional CaloCluster / CaloClusterMC pair:
// CaloClusterTruthMatch emits exactly one CaloClusterMC per CaloCluster
// in input order, so the two collections zip by index, and every
// CaloHitMC referenced by a CaloClusterMC inherits the corresponding
// CaloCluster::diskID(). The per-cluster hit lists themselves are NOT
// positionally matched (the MC list keeps only truth-matched hits and is
// re-sorted by MC energy), so no hit-by-hit pairing is attempted; only
// the disk, a cluster-level property, is taken. That is only sound for
// disk-local clusters, which is verified against the reco hits' crystals:
// a mixed-disk cluster (possible under cluster-association strategy 2,
// which does not enforce same-disk membership) throws rather than
// receiving a guessed disk. Configure the cluster tags only when the
// fallback is actually needed. CaloHitMC entries whose disk cannot be
// determined either way get null (unresolved) entrants.
//
// Grouping hits by entrant collapses secondary shower products
// (bremsstrahlung photons etc.) into their parent shower, recovering
// true shower membership for clustering truth definitions. Purity cuts
// and ambiguity handling are left to consumers.
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include "Offline/MCDataProducts/inc/CaloClusterMC.hh"
#include "Offline/MCDataProducts/inc/CaloHitEntrant.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"

#include <map>
#include <set>
#include <memory>
#include <utility>
#include <vector>

namespace mu2e {

  class CaloEntrantTruthMaker : public art::EDProducer
  {
     public:
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<art::InputTag> caloHitMCTag     { Name("caloHitMCTag"),     Comment("CaloHitMCCollection input tag") };
            fhicl::Atom<art::InputTag> caloClusterTag   { Name("caloClusterTag"),   Comment("CaloClusterCollection input tag (legacy disk resolution; empty disables)"), art::InputTag() };
            fhicl::Atom<art::InputTag> caloClusterMCTag { Name("caloClusterMCTag"), Comment("CaloClusterMCCollection input tag (legacy disk resolution; empty disables)"), art::InputTag() };
            fhicl::Atom<int>           diagLevel        { Name("diagLevel"),        Comment("Diagnostic level"), 0 };
        };
        using Parameters = art::EDProducer::Table<Config>;

        // member initialization follows declaration order below:
        // tokens, useClusterFallback_, diagLevel_
        explicit CaloEntrantTruthMaker(const Parameters& config) :
          EDProducer{config},
          caloHitMCToken_    {consumes<CaloHitMCCollection>(config().caloHitMCTag())},
          caloClusterToken_  {mayConsume<CaloClusterCollection>(config().caloClusterTag())},
          caloClusterMCToken_{mayConsume<CaloClusterMCCollection>(config().caloClusterMCTag())},
          useClusterFallback_(!config().caloClusterTag().empty()),
          diagLevel_         (config().diagLevel())
        {
            if (config().caloClusterTag().empty() != config().caloClusterMCTag().empty())
            {
                throw cet::exception("CALOENTRANT")
                    << "caloClusterTag and caloClusterMCTag must be both set or both empty\n";
            }
            produces<CaloHitEntrantCollection>();
        }

        void produce(art::Event& event) override;

     private:
        art::ProductToken<CaloHitMCCollection>     caloHitMCToken_;
        art::ProductToken<CaloClusterCollection>   caloClusterToken_;
        art::ProductToken<CaloClusterMCCollection> caloClusterMCToken_;
        bool useClusterFallback_;
        int  diagLevel_;
        bool warnedNoFallback_ = false;
  };

  //--------------------------------------------------------------------
  void CaloEntrantTruthMaker::produce(art::Event& event)
  {
      const auto  hitMCHandle = event.getValidHandle(caloHitMCToken_);
      const auto& hitMCs      = *hitMCHandle;
      const Calorimeter& cal  = *(GeomHandle<Calorimeter>());

      // Legacy disk resolution: every CaloHitMC referenced by a
      // CaloClusterMC inherits the disk of the corresponding CaloCluster.
      // Only the collection-level zip is an invariant of
      // CaloClusterTruthMatch (one CaloClusterMC per CaloCluster, input
      // order); the per-cluster hit lists are filtered and re-sorted, so
      // they are never paired positionally here.
      std::map<size_t, std::pair<int, size_t>> hitmcDisk; // HMC key -> (disk, first claiming cluster)
      if (useClusterFallback_)
      {
          const auto& clusters   = *event.getValidHandle(caloClusterToken_);
          const auto& clusterMCs = *event.getValidHandle(caloClusterMCToken_);
          if (clusters.size() != clusterMCs.size())
          {
              throw cet::exception("CALOENTRANT")
                  << "CaloCluster (" << clusters.size() << ") and CaloClusterMC ("
                  << clusterMCs.size() << ") collections are not parallel\n";
          }
          for (size_t ic = 0; ic < clusters.size(); ++ic)
          {
              const int disk = clusters[ic].diskID();

              // The cluster-wide disk is only meaningful for a disk-local
              // cluster. Association strategy 2 (ClusterAssociator) can merge
              // proto-clusters across disks while CaloCluster keeps only the
              // seed's disk; the filtered, re-sorted MC list cannot recover
              // per-hit disks for such a cluster, so fail closed instead of
              // assigning a guessed disk.
              for (const auto& hit : clusters[ic].caloHitsPtrVector())
              {
                  if (hit.isNull()) continue;
                  const int hitDisk = cal.crystal(hit->crystalID()).diskID();
                  if (hitDisk != disk)
                  {
                      throw cet::exception("CALOENTRANT")
                          << "cluster " << ic << " is not disk-local: diskID " << disk
                          << " but hit crystal " << hit->crystalID() << " is on disk "
                          << hitDisk << "; the cluster-wide disk fallback cannot label it\n";
                  }
              }

              for (const auto& hmc : clusterMCs[ic].caloHitMCs())
              {
                  if (hmc.isNull()) continue;
                  if (hmc.id() != hitMCHandle.id())
                  {
                      throw cet::exception("CALOENTRANT")
                          << "CaloClusterMC references CaloHitMCCollection " << hmc.id()
                          << " but caloHitMCTag resolves to " << hitMCHandle.id()
                          << "; caloClusterMCTag and caloHitMCTag are inconsistent\n";
                  }
                  if (hmc.key() >= hitMCs.size())
                  {
                      throw cet::exception("CALOENTRANT")
                          << "CaloClusterMC " << ic << " references CaloHitMC key " << hmc.key()
                          << " outside collection size " << hitMCs.size() << "\n";
                  }
                  const auto ins = hitmcDisk.emplace(hmc.key(), std::make_pair(disk, ic));
                  if (!ins.second && ins.first->second.first != disk)
                  {
                      throw cet::exception("CALOENTRANT")
                          << "CaloHitMC key " << hmc.key() << " (ProductID " << hmc.id()
                          << ") is claimed by cluster " << ins.first->second.second
                          << " on disk " << ins.first->second.first
                          << " and by cluster " << ic << " on disk " << disk << "\n";
                  }
              }
          }
      }

      // Disk of a CaloHitMC: from the crystalID member when filled (newer
      // productions), else the cluster lookup (legacy files), else -1.
      auto diskOf = [&](size_t idx, const CaloHitMC& chmc) -> int
      {
          if (chmc.crystalID() >= 0) return cal.crystal(chmc.crystalID()).diskID();
          const auto found = hitmcDisk.find(idx);
          return (found != hitmcDisk.end()) ? found->second.first : -1;
      };

      // Pass 1: record which disks each SimParticle deposited in. Keyed by
      // art::Ptr (ProductID + key) so deposits referencing more than one
      // SimParticle collection cannot collide.
      std::map<art::Ptr<SimParticle>, std::set<int>> simDisks;
      for (size_t idx = 0; idx < hitMCs.size(); ++idx)
      {
          const int disk = diskOf(idx, hitMCs[idx]);
          if (disk < 0) continue;
          for (const auto& edep : hitMCs[idx].energyDeposits())
          {
              if (edep.sim().isNonnull()) simDisks[edep.sim()].insert(disk);
          }
      }

      // Pass 2: for each deposit, walk the parent chain upward; the
      // entrant is the highest ancestor that also deposited in the same
      // disk. The walk stops at the chain root; a parent Ptr that cannot
      // be read back (ancestor collection dropped from the file) throws.
      auto output = std::make_unique<CaloHitEntrantCollection>();
      output->reserve(hitMCs.size());
      std::map<std::pair<art::Ptr<SimParticle>, int>, art::Ptr<SimParticle>> entrantCache;
      size_t nUnresolved(0);

      for (size_t idx = 0; idx < hitMCs.size(); ++idx)
      {
          const auto& chmc = hitMCs[idx];
          std::vector<art::Ptr<SimParticle>> entrants;

          const int disk = diskOf(idx, chmc);
          if (disk < 0)
          {
              // no crystalID member and not referenced by any cluster:
              // entrants unresolvable; keep the collection index-parallel
              // with null entries
              entrants.assign(chmc.nParticles(), art::Ptr<SimParticle>());
              ++nUnresolved;
              output->emplace_back(art::Ptr<CaloHitMC>(hitMCHandle, idx), std::move(entrants));
              continue;
          }
          entrants.reserve(chmc.nParticles());

          for (const auto& edep : chmc.energyDeposits())
          {
              const art::Ptr<SimParticle>& sim = edep.sim();
              if (sim.isNull()) { entrants.emplace_back(); continue; }

              const auto key = std::make_pair(sim, disk);
              auto it = entrantCache.find(key);
              if (it == entrantCache.end())
              {
                  if (!sim.isAvailable())
                  {
                      throw cet::exception("CALOENTRANT")
                          << "CaloEntrantTruthMaker: depositor SimParticle (ProductID " << sim.id()
                          << ", key " << sim.key() << ") cannot be read from this file"
                          << " - its SimParticle collection was dropped;"
                          << " cannot walk the parent chain\n";
                  }
                  art::Ptr<SimParticle> entrant = sim;
                  for (auto p = sim->parent(); p.isNonnull(); p = p->parent())
                  {
                      if (!p.isAvailable())
                      {
                          throw cet::exception("CALOENTRANT")
                              << "CaloEntrantTruthMaker: SimParticle parent (ProductID " << p.id()
                              << ", key " << p.key() << ") cannot be read from this file"
                              << " - an ancestor SimParticle collection was dropped;"
                              << " cannot walk the parent chain\n";
                      }
                      const auto fnd = simDisks.find(p);
                      if (fnd != simDisks.end() && fnd->second.count(disk)) entrant = p;
                  }
                  it = entrantCache.emplace(key, entrant).first;
              }
              entrants.push_back(it->second);
          }
          output->emplace_back(art::Ptr<CaloHitMC>(hitMCHandle, idx), std::move(entrants));
      }

      if (nUnresolved > 0 && !useClusterFallback_ && !warnedNoFallback_)
      {
          warnedNoFallback_ = true;
          mf::LogWarning("CaloEntrantTruthMaker")
              << nUnresolved << " of " << hitMCs.size() << " CaloHitMC have no crystalID"
              << " and no cluster fallback is configured; their entrants are null."
              << " Set caloClusterTag/caloClusterMCTag to resolve disks on legacy files.";
      }
      if (diagLevel_ > 0)
      {
          mf::LogInfo("CaloEntrantTruthMaker") << "produced " << output->size()
                                               << " CaloHitEntrant for " << hitMCs.size() << " CaloHitMC"
                                               << " (" << hitmcDisk.size() << " cluster-referenced, "
                                               << nUnresolved << " unresolved)";
      }

      event.put(std::move(output));
  }
}

DEFINE_ART_MODULE(mu2e::CaloEntrantTruthMaker)
