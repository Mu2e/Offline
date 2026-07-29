//
// Assign each calorimeter MC energy deposit to its calo-entrant ancestor:
// the highest SimParticle in the Geant4 parent chain that also deposited
// energy in the same calorimeter disk (the shower originator).
//
// For each CaloHitMC in the input collection this module emits a
// CaloHitEntrant whose entrants vector is aligned with
// CaloHitMC::energyDeposits(); the output collection is index-parallel
// to the input CaloHitMCCollection.
//
// Crystal (hence disk) resolution prefers CaloHitMC::crystalID() when it
// is filled. For files produced before that member existed (reads back
// -1) the crystal is resolved through the CaloCluster <-> CaloClusterMC
// pairing (index-parallel collections with positionally-matched hit
// lists, the invariant CaloClusterTruthMatch establishes). CaloHitMC
// entries whose crystal cannot be determined either way get null
// (unresolved) entrants.
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

namespace mu2e {

  class CaloEntrantTruthMaker : public art::EDProducer
  {
     public:
        struct Config
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;
            fhicl::Atom<art::InputTag> caloHitMCTag     { Name("caloHitMCTag"),     Comment("CaloHitMCCollection input tag") };
            fhicl::Atom<art::InputTag> caloClusterTag   { Name("caloClusterTag"),   Comment("CaloClusterCollection input tag (legacy crystal resolution)") };
            fhicl::Atom<art::InputTag> caloClusterMCTag { Name("caloClusterMCTag"), Comment("CaloClusterMCCollection input tag (legacy crystal resolution)") };
            fhicl::Atom<int>           diagLevel        { Name("diagLevel"),        Comment("Diagnostic level"), 0 };
        };

        explicit CaloEntrantTruthMaker(const art::EDProducer::Table<Config>& config) :
          EDProducer{config},
          caloHitMCToken_    {consumes<CaloHitMCCollection>(config().caloHitMCTag())},
          caloClusterToken_  {consumes<CaloClusterCollection>(config().caloClusterTag())},
          caloClusterMCToken_{consumes<CaloClusterMCCollection>(config().caloClusterMCTag())},
          diagLevel_         (config().diagLevel())
        {
            produces<CaloHitEntrantCollection>();
        }

        void produce(art::Event& event) override;

     private:
        art::ProductToken<CaloHitMCCollection>     caloHitMCToken_;
        art::ProductToken<CaloClusterCollection>   caloClusterToken_;
        art::ProductToken<CaloClusterMCCollection> caloClusterMCToken_;
        int diagLevel_;
  };

  //--------------------------------------------------------------------
  void CaloEntrantTruthMaker::produce(art::Event& event)
  {
      const auto& hitMCs     = *event.getValidHandle(caloHitMCToken_);
      const auto& clusters   = *event.getValidHandle(caloClusterToken_);
      const auto& clusterMCs = *event.getValidHandle(caloClusterMCToken_);
      const Calorimeter& cal = *(GeomHandle<Calorimeter>());

      // Legacy crystal resolution: CaloCluster <-> CaloClusterMC pairing.
      // CaloClusterTruthMatch emits exactly one CaloClusterMC per
      // CaloCluster in input order with one CaloHitMC per hit, so any
      // cardinality mismatch is a wiring error, not a soft condition.
      if (clusters.size() != clusterMCs.size())
      {
          throw cet::exception("CALOENTRANT")
              << "CaloCluster (" << clusters.size() << ") and CaloClusterMC ("
              << clusterMCs.size() << ") collections are not parallel\n";
      }
      std::map<size_t, int> hitmcCrystal;
      for (size_t ic = 0; ic < clusters.size(); ++ic)
      {
          const auto& hits   = clusters[ic].caloHitsPtrVector();
          const auto& hitmcs = clusterMCs[ic].caloHitMCs();
          if (hits.size() != hitmcs.size())
          {
              throw cet::exception("CALOENTRANT")
                  << "cluster " << ic << ": " << hits.size() << " CaloHits vs "
                  << hitmcs.size() << " CaloHitMCs\n";
          }
          for (size_t ih = 0; ih < hits.size(); ++ih)
          {
              if (hits[ih].isNull() || hitmcs[ih].isNull())
              {
                  throw cet::exception("CALOENTRANT")
                      << "cluster " << ic << " hit " << ih << ": null Ptr in cluster hit lists\n";
              }
              hitmcCrystal[hitmcs[ih].key()] = hits[ih]->crystalID();
          }
      }

      // Crystal of a CaloHitMC: the member when filled (newer
      // productions), else the cluster-pairing lookup (legacy files).
      auto crystalOf = [&](size_t idx, const CaloHitMC& chmc) -> int
      {
          if (chmc.crystalID() >= 0) return chmc.crystalID();
          const auto found = hitmcCrystal.find(idx);
          return (found != hitmcCrystal.end()) ? found->second : -1;
      };

      // Pass 1: record which disks each SimParticle deposited in.
      // Keyed by SimParticle id to match the id convention exposed to
      // ntuple consumers (SimParticle::id().asInt()).
      std::map<int, std::set<int>> simDisks;
      for (size_t idx = 0; idx < hitMCs.size(); ++idx)
      {
          const int crystal = crystalOf(idx, hitMCs[idx]);
          if (crystal < 0) continue;
          int disk = cal.crystal(crystal).diskID();
          for (const auto& edep : hitMCs[idx].energyDeposits())
          {
              if (edep.sim().isNonnull()) simDisks[edep.sim()->id().asInt()].insert(disk);
          }
      }

      // Pass 2: for each deposit, walk the parent chain upward; the
      // entrant is the highest ancestor that also deposited in the same
      // disk. The walk stops at the chain root or at an unresolvable
      // parent pointer (event compression).
      auto output = std::make_unique<CaloHitEntrantCollection>();
      output->reserve(hitMCs.size());
      std::map<std::pair<int,int>, art::Ptr<SimParticle>> entrantCache; // (simId, disk) -> entrant

      for (size_t idx = 0; idx < hitMCs.size(); ++idx)
      {
          const auto& chmc = hitMCs[idx];
          CaloHitEntrant che;

          const int crystal = crystalOf(idx, chmc);
          if (crystal < 0)
          {
              // not referenced by any cluster and no crystal member:
              // entrants unresolvable; keep the collection index-parallel
              // with null entries
              che.entrants.assign(chmc.nParticles(), art::Ptr<SimParticle>());
              output->push_back(std::move(che));
              continue;
          }
          const int disk = cal.crystal(crystal).diskID();
          che.entrants.reserve(chmc.nParticles());

          for (const auto& edep : chmc.energyDeposits())
          {
              const art::Ptr<SimParticle>& sim = edep.sim();
              if (sim.isNull()) { che.entrants.emplace_back(); continue; }

              const auto key = std::make_pair(sim->id().asInt(), disk);
              auto it = entrantCache.find(key);
              if (it == entrantCache.end())
              {
                  art::Ptr<SimParticle> entrant = sim;
                  for (auto p = sim->parent(); p.isNonnull(); p = p->parent())
                  {
                      auto fnd = simDisks.find(p->id().asInt());
                      if (fnd != simDisks.end() && fnd->second.count(disk)) entrant = p;
                  }
                  it = entrantCache.emplace(key, entrant).first;
              }
              che.entrants.push_back(it->second);
          }
          output->push_back(std::move(che));
      }

      if (diagLevel_ > 0)
      {
          mf::LogInfo("CaloEntrantTruthMaker") << "produced " << output->size()
                                               << " CaloHitEntrant for " << hitMCs.size() << " CaloHitMC"
                                               << " (" << hitmcCrystal.size() << " cluster-referenced)";
      }

      event.put(std::move(output));
  }
}

DEFINE_ART_MODULE(mu2e::CaloEntrantTruthMaker)
