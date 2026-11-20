//
// CaloHitGraphParityDump — dumps per-event CaloHit metadata and the
// GNN cluster assignments to a flat TTree, so a Python comparison
// script can replay the same events through the Python pipeline and
// assert agreement on the graph and the cluster labels.
//
// Output TTree (in the TFileService output), one entry per event-disk:
//
//   eventID              : ULong64_t
//   diskID               : Int_t
//   nHits                : Int_t                     (CaloHits on this disk)
//   crystalID[nHits]     : std::vector<int>          (per hit)
//   time_ns[nHits]       : std::vector<float>
//   eDep_MeV[nHits]      : std::vector<float>
//   xLocal_mm[nHits]     : std::vector<double>       (crystal position in the
//   yLocal_mm[nHits]     : std::vector<double>        disk frame, as the graph
//                                                     maker uses it)
//   gnnLabel[nHits]      : std::vector<int>          (-1 if the hit was dropped
//                                                     by the minHits / minEnergyMeV
//                                                     cleanup; 0..K-1 otherwise)
//   gnnHitRank[nHits]    : std::vector<int>          (position of the hit in its
//                                                     cluster's hit list; -1 if
//                                                     dropped)
//
// With caloHitGraphCollection set, also the graph the cluster maker saw,
// with node indices in the order of the hits above:
//
//   edgeSrc[nEdges], edgeDst[nEdges] : std::vector<int>
//   nodeX[6 * nHits]                 : std::vector<float>  (normalized)
//
// Used for the C++/Python parity check of the graphs and the clusters.
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/DataProducts/inc/CrystalId.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloHitGraph.hh"

#include "TTree.h"

#include <cstdint>
#include <map>
#include <optional>
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
        fhicl::OptionalAtom<art::InputTag> caloHitGraphCollection {
          Name("caloHitGraphCollection"),
          Comment("Optional CaloHitGraph collection; when set, its edges and node features are dumped") };
      };

      explicit CaloHitGraphParityDump(const art::EDAnalyzer::Table<Config>& config);

      void analyze(const art::Event& event) override;
      void beginJob() override;

    private:
      art::ProductToken<CaloHitCollection>                      hitToken_;
      art::ProductToken<CaloClusterCollection>                  clusterToken_;
      std::optional<art::ProductToken<CaloHitGraphCollection>>  graphToken_;

      TTree*                  tree_  = nullptr;
      std::uint64_t           bEventID_  = 0;
      int                     bDiskID_   = -1;
      int                     bNHits_    = 0;
      std::vector<int>        bCrystalID_;
      std::vector<float>      bTime_;
      std::vector<float>      bEDep_;
      std::vector<double>     bXLocal_;
      std::vector<double>     bYLocal_;
      std::vector<int>        bGnnLabel_;
      std::vector<int>        bGnnHitRank_;
      std::vector<int>        bEdgeSrc_;
      std::vector<int>        bEdgeDst_;
      std::vector<float>      bNodeX_;
  };


  CaloHitGraphParityDump::CaloHitGraphParityDump(
      const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config},
    hitToken_    {consumes<CaloHitCollection>(config().caloHitCollection())},
    clusterToken_{consumes<CaloClusterCollection>(config().caloClusterCollection())}
  {
    art::InputTag graphTag;
    if (config().caloHitGraphCollection(graphTag)) {
      graphToken_ = consumes<CaloHitGraphCollection>(graphTag);
    }
  }


  void CaloHitGraphParityDump::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    tree_ = tfs->make<TTree>("parity", "GNN-cluster parity dump");
    tree_->Branch("eventID",    &bEventID_);
    tree_->Branch("diskID",     &bDiskID_);
    tree_->Branch("nHits",      &bNHits_);
    tree_->Branch("crystalID",  &bCrystalID_);
    tree_->Branch("time_ns",    &bTime_);
    tree_->Branch("eDep_MeV",   &bEDep_);
    tree_->Branch("xLocal_mm",  &bXLocal_);
    tree_->Branch("yLocal_mm",  &bYLocal_);
    tree_->Branch("gnnLabel",   &bGnnLabel_);
    tree_->Branch("gnnHitRank", &bGnnHitRank_);
    if (graphToken_) {
      tree_->Branch("edgeSrc",  &bEdgeSrc_);
      tree_->Branch("edgeDst",  &bEdgeDst_);
      tree_->Branch("nodeX",    &bNodeX_);
    }
  }


  void CaloHitGraphParityDump::analyze(const art::Event& event)
  {
    auto hitsHandle    = event.getHandle<CaloHitCollection>(hitToken_);
    auto clustersHandle= event.getHandle<CaloClusterCollection>(clusterToken_);
    if (!hitsHandle.isValid() || !clustersHandle.isValid()) return;

    const Calorimeter& cal = *(GeomHandle<Calorimeter>());
    const auto& hits     = *hitsHandle;
    const auto& clusters = *clustersHandle;

    // Per-CaloHit-pointer cluster id and rank. The cluster id is its
    // position within the per-disk subset of the CaloClusterCollection
    // (each disk-graph is reconstructed independently); the rank is the
    // hit's position in that cluster's hit list.
    std::map<int, std::map<const CaloHit*, std::pair<int, int>>> diskClusterByHit;
    std::map<int, int> nextLabelByDisk;
    for (const auto& cluster : clusters) {
      const int disk = cluster.diskID();
      const int label = nextLabelByDisk[disk]++;
      int rank = 0;
      for (const auto& hp : cluster.caloHitsPtrVector()) {
        diskClusterByHit[disk][hp.get()] = {label, rank++};
      }
    }

    // Graphs by disk, if requested.
    std::map<int, const CaloHitGraph*> graphByDisk;
    if (graphToken_) {
      auto graphHandle = event.getHandle<CaloHitGraphCollection>(*graphToken_);
      if (!graphHandle.isValid()) return;
      for (const auto& g : *graphHandle) graphByDisk[g.diskID] = &g;
    }

    // Partition CaloHits by disk and emit one TTree entry per disk.
    std::map<int, std::vector<const CaloHit*>> hitsByDisk;
    for (const auto& h : hits) {
      if (!CrystalId(h.crystalID()).isValid()) continue; // pin diodes never enter a graph
      const int d = cal.crystal(h.crystalID()).diskID();
      hitsByDisk[d].push_back(&h);
    }

    bEventID_ = static_cast<std::uint64_t>(event.id().event());
    for (const auto& [disk, diskHits] : hitsByDisk) {
      bDiskID_   = disk;
      bNHits_    = static_cast<int>(diskHits.size());
      bCrystalID_.clear();  bTime_.clear();  bEDep_.clear();
      bXLocal_.clear();     bYLocal_.clear();
      bGnnLabel_.clear();   bGnnHitRank_.clear();
      bEdgeSrc_.clear();    bEdgeDst_.clear();  bNodeX_.clear();
      const auto& hitMap = diskClusterByHit[disk];
      std::map<const CaloHit*, int> indexOf;
      for (const CaloHit* h : diskHits) {
        indexOf[h] = static_cast<int>(bCrystalID_.size());
        const auto& pos = cal.crystal(h->crystalID()).localPosition();
        bCrystalID_.push_back(h->crystalID());
        bTime_     .push_back(h->time());
        bEDep_     .push_back(h->energyDep());
        bXLocal_   .push_back(pos.x());
        bYLocal_   .push_back(pos.y());
        auto it = hitMap.find(h);
        bGnnLabel_  .push_back(it != hitMap.end() ? it->second.first  : -1);
        bGnnHitRank_.push_back(it != hitMap.end() ? it->second.second : -1);
      }

      if (graphToken_) {
        auto git = graphByDisk.find(disk);
        if (git == graphByDisk.end()) {
          throw cet::exception("CaloHitGraphParityDump")
            << "no CaloHitGraph for disk " << disk << " in event " << bEventID_;
        }
        const CaloHitGraph& g = *git->second;
        if (g.nNodes != bNHits_) {
          throw cet::exception("CaloHitGraphParityDump")
            << "graph of disk " << disk << " has " << g.nNodes
            << " nodes, the disk has " << bNHits_ << " hits";
        }
        // Graph node i -> dump hit index, through the CaloHit pointers.
        std::vector<int> dumpIndex(g.nNodes);
        bNodeX_.assign(6 * g.nNodes, 0.0f);
        for (int i = 0; i < g.nNodes; ++i) {
          auto ix = indexOf.find(g.caloHitPtrs[i].get());
          if (ix == indexOf.end()) {
            throw cet::exception("CaloHitGraphParityDump")
              << "graph node " << i << " of disk " << disk << " is not a dumped hit";
          }
          dumpIndex[i] = ix->second;
          for (int k = 0; k < 6; ++k) bNodeX_[6 * ix->second + k] = g.x[6 * i + k];
        }
        for (int e = 0; e < g.nEdges; ++e) {
          bEdgeSrc_.push_back(dumpIndex[g.edgeIndex[e]]);
          bEdgeDst_.push_back(dumpIndex[g.edgeIndex[g.nEdges + e]]);
        }
      }
      tree_->Fill();
    }
  }

}

DEFINE_ART_MODULE(mu2e::CaloHitGraphParityDump)
