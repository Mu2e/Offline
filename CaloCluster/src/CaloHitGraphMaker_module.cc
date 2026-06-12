//
// CaloHitGraphMaker — first half of the GNN clustering split design
// (see calorimeter/GNN/docs/offline_integration.md §1.2).
//
// Per event, partitions the input CaloHitCollection by disk and runs
// the C++ port of the Python graph builder (GnnGraphBuilder) once per
// disk. Emits a CaloHitGraphCollection with one entry per disk that
// has at least one CaloHit. The emitted graphs already carry
// z-score-normalised feature tensors so the downstream cluster module
// has no normalisation responsibility.
//
// Frozen hyper-parameters and feature column order live in the
// canonical norm sidecar (calorimeter/GNN/scripts/export_norm_stats.py)
// and are validated against canonical names at sidecar load time.
// FHiCL parameters expose the graph-construction knobs so swapping
// future models with different feature lists is config-driven.
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/CaloCluster/inc/GnnGraphBuilder.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloHitGraph.hh"

#include <memory>
#include <string>
#include <vector>


namespace mu2e {

  class CaloHitGraphMaker : public art::EDProducer
  {
    public:
      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> caloHitCollection {
          Name("caloHitCollection"),
          Comment("CaloHit collection to read") };
        fhicl::Atom<std::string>   normSidecar {
          Name("normSidecar"),
          Comment("Relative path to the JSON norm sidecar; resolved by ConfigFileLookupPolicy") };
        fhicl::Atom<double>        rMax {
          Name("rMax"),
          Comment("Spatial radius cut for the radius graph (mm)"),
          210.0 };
        fhicl::Atom<double>        dtMax {
          Name("dtMax"),
          Comment("Maximum |dt| between connected hits (ns)"),
          25.0 };
        fhicl::Atom<unsigned>      kMin {
          Name("kMin"),
          Comment("kNN-fallback floor for nodes with degree < kMin"),
          3 };
        fhicl::Atom<unsigned>      kMax {
          Name("kMax"),
          Comment("Per-source-node degree cap"),
          20 };
      };

      explicit CaloHitGraphMaker(const art::EDProducer::Table<Config>& config) :
        art::EDProducer{config},
        caloHitToken_{consumes<CaloHitCollection>(config().caloHitCollection())},
        builder_(makeBuilder(config()))
      {
        produces<CaloHitGraphCollection>();
      }

      void produce(art::Event& event) override;

    private:
      static GnnGraphBuilder makeBuilder(const Config& cfg);

      art::ProductToken<CaloHitCollection> caloHitToken_;
      GnnGraphBuilder                      builder_;
  };


  GnnGraphBuilder
  CaloHitGraphMaker::makeBuilder(const Config& cfg)
  {
    GnnGraphBuilder::Config gcfg;
    gcfg.rMax  = cfg.rMax();
    gcfg.dtMax = cfg.dtMax();
    gcfg.kMin  = cfg.kMin();
    gcfg.kMax  = cfg.kMax();
    ConfigFileLookupPolicy lookup;
    const std::string statsPath = lookup(cfg.normSidecar());
    auto stats = GnnGraphBuilder::loadStatsFromJson(statsPath);
    return GnnGraphBuilder(gcfg, stats);
  }


  void CaloHitGraphMaker::produce(art::Event& event)
  {
    auto hitsHandle = event.getHandle<CaloHitCollection>(caloHitToken_);
    auto out = std::make_unique<CaloHitGraphCollection>();

    const auto& hits = *hitsHandle;
    if (hits.empty()) {
      event.put(std::move(out));
      return;
    }

    const Calorimeter& cal = *(GeomHandle<Calorimeter>());

    // Partition CaloHits by disk. Two passes so nDisks isn't hardcoded:
    // first pass to discover the disk set, then bucket pointers + Ptrs.
    std::vector<int> diskIDs;
    std::vector<std::vector<const CaloHit*>>      hitsByDisk;
    std::vector<std::vector<art::Ptr<CaloHit>>>   ptrsByDisk;

    for (std::size_t i = 0; i < hits.size(); ++i) {
      const CaloHit& h = hits[i];
      const int disk = cal.crystal(h.crystalID()).diskID();

      // Find/register the bucket for this disk.
      std::size_t bucket = diskIDs.size();
      for (std::size_t k = 0; k < diskIDs.size(); ++k) {
        if (diskIDs[k] == disk) { bucket = k; break; }
      }
      if (bucket == diskIDs.size()) {
        diskIDs.push_back(disk);
        hitsByDisk.emplace_back();
        ptrsByDisk.emplace_back();
      }
      hitsByDisk[bucket].push_back(&h);
      ptrsByDisk[bucket].emplace_back(hitsHandle, i);
    }

    // Build one graph per disk that has hits.
    out->resize(diskIDs.size());
    for (std::size_t k = 0; k < diskIDs.size(); ++k) {
      builder_.buildGraph(diskIDs[k],
                          hitsByDisk[k], ptrsByDisk[k],
                          cal, (*out)[k]);
    }

    event.put(std::move(out));
  }

}

DEFINE_ART_MODULE(mu2e::CaloHitGraphMaker)
