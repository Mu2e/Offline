#ifndef CaloCluster_GnnGraphBuilder_hh
#define CaloCluster_GnnGraphBuilder_hh
//
// C++ port of calorimeter/GNN/src/data/graph_builder.py.
//
// Builds one CaloHitGraph per calorimeter disk per event:
//   1. Collect CaloHits per disk; look up (x, y) in the disk-local
//      frame from the Calorimeter geometry service.
//   2. Brute-force pairwise distance loop with r_max cut for the
//      radius graph (faithful to scipy.spatial.cKDTree.query_pairs).
//   3. Time filter |dt| < dt_max ns.
//   4. kNN fallback for nodes with degree < k_min after the radius+time pass.
//   5. Per-source-node degree cap at k_max (keep the k_max nearest dsts).
//   6. Compute 6 node features and 8 edge features.
//   7. Z-score normalise using the train-split stats from the JSON
//      sidecar passed at construction (loaded via loadStatsFromJson).
//
// Feature column order is canonical and matches the model's
// metadata_props (see calorimeter/GNN/docs/onnx_deployment.md):
//
//   nodes : log_e, t, x, y, r, e_rel
//   edges : dx, dy, d, dt, dlog_e, asym_e, logsum_e, dr
//

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloHitGraph.hh"

#include "canvas/Persistency/Common/Ptr.h"

#include <string>
#include <vector>

namespace mu2e {

  class GnnGraphBuilder
  {
    public:
      // Per-feature normalisation stats (z-score: (x - mean) / std).
      struct Stats
      {
        std::vector<float> nodeMean;  // size 6
        std::vector<float> nodeStd;   // size 6
        std::vector<float> edgeMean;  // size 8
        std::vector<float> edgeStd;   // size 8
      };

      struct Config
      {
        double rMax   = 210.0;  // mm — radius graph cut
        double dtMax  = 25.0;   // ns — time-coincidence cut
        unsigned kMin = 3;      // kNN fallback floor
        unsigned kMax = 20;     // per-source-node degree cap
      };

      GnnGraphBuilder(const Config& cfg, const Stats& stats)
        : cfg_(cfg), stats_(stats) {}

      // Load Stats from the JSON sidecar produced by
      // calorimeter/GNN/scripts/export_norm_stats.py. Throws
      // cet::exception on missing keys, wrong sizes, or canonical
      // node/edge feature-name mismatches.
      static Stats loadStatsFromJson(const std::string& jsonPath);

      // Build one CaloHitGraph for one disk.
      //   diskID    — destination disk for the emitted graph
      //   hits      — pointers to the CaloHits on this disk
      //   ptrs      — art::Ptr back to each hit, parallel to `hits`
      //   cal       — geometry handle for crystal positions
      //   out       — populated in place (cleared first)
      void buildGraph(int diskID,
                      const std::vector<const CaloHit*>& hits,
                      const std::vector<art::Ptr<CaloHit>>& ptrs,
                      const Calorimeter& cal,
                      CaloHitGraph& out) const;

    private:
      Config cfg_;
      Stats  stats_;
  };

}

#endif
