//
// CaloHitGraphParityDump — dumps per-event CaloHit metadata and the
// GNN cluster assignments to a flat TTree, so a Python comparison
// script can replay the same events through the Python pipeline and
// assert byte-exact agreement on cluster labels.
//
// Output TTree (in the TFileService output):
//
//   eventID              : ULong64_t
//   diskID               : Int_t
//   nHits                : Int_t                     (CaloHits on this disk)
//   crystalID[nHits]     : std::vector<int>          (per hit)
//   time_ns[nHits]       : std::vector<float>
//   eDep_MeV[nHits]      : std::vector<float>
//   gnnLabel[nHits]      : std::vector<int>          (-1 if hit was dropped
//                                                     by min_hits / min_E cut;
//                                                     0..K-1 otherwise)
//
// One entry per event-disk. The Python comparison script reads the
// TTree, rebuilds the per-disk graph from the dumped CaloHit info,
// runs the Python CaloClusterNet + cluster_reco, and asserts that
// label vectors match.
//
// Used as the Stage 3 parity check in Task 16g.
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"

#include "TTree.h"

#include <cstdint>
#include <map>
#include <vector>


namespace mu2e {

  class CaloHitGraphParityDump : public art::EDAnalyzer
  {
    public:
      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> caloHitCollection {
          Name("caloHitCollection"),
          Comment("CaloHit collection (input to the GNN graph maker)") };
        fhicl::Atom<art::InputTag> caloClusterCollection {
          Name("caloClusterCollection"),
          Comment("CaloCluster collection emitted by CaloClusterMakerGNN, e.g. caloClusterMakerGNN:GNN") };
      };

      explicit CaloHitGraphParityDump(const art::EDAnalyzer::Table<Config>& config);

      void analyze(const art::Event& event) override;
      void beginJob() override;

    private:
      art::ProductToken<CaloHitCollection>       hitToken_;
      art::ProductToken<CaloClusterCollection>   clusterToken_;

      TTree*                  tree_  = nullptr;
      std::uint64_t           bEventID_  = 0;
      int                     bDiskID_   = -1;
      int                     bNHits_    = 0;
      std::vector<int>        bCrystalID_;
      std::vector<float>      bTime_;
      std::vector<float>      bEDep_;
      std::vector<int>        bGnnLabel_;
  };


  CaloHitGraphParityDump::CaloHitGraphParityDump(
      const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    hitToken_    {consumes<CaloHitCollection>(config().caloHitCollection())},
    clusterToken_{consumes<CaloClusterCollection>(config().caloClusterCollection())}
  {
  }


  void CaloHitGraphParityDump::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    tree_ = tfs->make<TTree>("parity", "GNN-cluster parity dump");
    tree_->Branch("eventID",   &bEventID_);
    tree_->Branch("diskID",    &bDiskID_);
    tree_->Branch("nHits",     &bNHits_);
    tree_->Branch("crystalID", &bCrystalID_);
    tree_->Branch("time_ns",   &bTime_);
    tree_->Branch("eDep_MeV",  &bEDep_);
    tree_->Branch("gnnLabel",  &bGnnLabel_);
  }


  void CaloHitGraphParityDump::analyze(const art::Event& event)
  {
    auto hitsHandle    = event.getHandle<CaloHitCollection>(hitToken_);
    auto clustersHandle= event.getHandle<CaloClusterCollection>(clusterToken_);
    if (!hitsHandle.isValid() || !clustersHandle.isValid()) return;

    const Calorimeter& cal = *(GeomHandle<Calorimeter>());
    const auto& hits     = *hitsHandle;
    const auto& clusters = *clustersHandle;

    // Build per-CaloHit-pointer → cluster-id-on-this-disk map.
    // The cluster index is its position within the per-disk subset of
    // the CaloClusterCollection (matching what the Python pipeline
    // emits, since each disk-graph is reconstructed independently).
    std::map<int, std::map<const CaloHit*, int>> diskClusterByHit;  // disk -> hit* -> labelIdx
    std::map<int, int> nextLabelByDisk;
    for (const auto& cluster : clusters) {
      const int disk = cluster.diskID();
      const int label = nextLabelByDisk[disk]++;
      for (const auto& hp : cluster.caloHitsPtrVector()) {
        diskClusterByHit[disk][hp.get()] = label;
      }
    }

    // Partition CaloHits by disk and emit one TTree entry per disk.
    std::map<int, std::vector<const CaloHit*>> hitsByDisk;
    for (const auto& h : hits) {
      const int d = cal.crystal(h.crystalID()).diskID();
      hitsByDisk[d].push_back(&h);
    }

    bEventID_ = static_cast<std::uint64_t>(event.id().event());
    for (const auto& [disk, diskHits] : hitsByDisk) {
      bDiskID_   = disk;
      bNHits_    = static_cast<int>(diskHits.size());
      bCrystalID_.clear();  bCrystalID_.reserve(diskHits.size());
      bTime_    .clear();   bTime_.reserve(diskHits.size());
      bEDep_    .clear();   bEDep_.reserve(diskHits.size());
      bGnnLabel_.clear();   bGnnLabel_.reserve(diskHits.size());
      const auto& hitMap = diskClusterByHit[disk];
      for (const CaloHit* h : diskHits) {
        bCrystalID_.push_back(h->crystalID());
        bTime_    .push_back(h->time());
        bEDep_    .push_back(h->energyDep());
        auto it = hitMap.find(h);
        bGnnLabel_.push_back(it != hitMap.end() ? it->second : -1);
      }
      tree_->Fill();
    }
  }

}

DEFINE_ART_MODULE(mu2e::CaloHitGraphParityDump)
