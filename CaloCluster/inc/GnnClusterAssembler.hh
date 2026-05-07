#ifndef CaloCluster_GnnClusterAssembler_hh
#define CaloCluster_GnnClusterAssembler_hh
//
// C++ port of calorimeter/GNN/src/inference/cluster_reco.py for the
// CCN+BFS10 recipe (the winning configuration in
// docs/findings.md §7.4).
//
// Steps applied to the directed edge logits emitted by the ONNX model:
//   1. Sigmoid → per-edge probabilities.
//   2. Symmetrise: for each unordered pair {i, j}, take the mean of
//      p_ij and p_ji.
//   3. Threshold at tauEdge.
//   4. BFS traversal seeded from highest-energy hits — hits with
//      energy >= bfsExpandCut continue the BFS; lower-energy hits join
//      but cannot recruit. Mirrors Offline's ClusterFinder ExpandCut.
//   5. Cleanup: drop clusters with fewer than minHits hits or total
//      energy below minEnergyMeV.
//   6. Relabel to contiguous IDs.
//
// Returns labels[N] where labels[i] = cluster ID >= 0 or -1 (dropped).
//

#include <cstdint>
#include <vector>

namespace mu2e {

  class GnnClusterAssembler
  {
    public:
      struct Config
      {
        double   tauEdge       = 0.20;  // probability threshold (model-specific)
        double   bfsExpandCut  = 10.0;  // MeV — BFS-style ExpandCut
        unsigned minHits       = 2;     // drop clusters smaller than this
        double   minEnergyMeV  = 10.0;  // drop clusters below this total energy
      };

      explicit GnnClusterAssembler(const Config& cfg) : cfg_(cfg) {}

      // nNodes        : number of hits in the graph
      // edgeIndex     : flat (2 * E) int64s, src row first then dst row,
      //                 matching the CaloHitGraph layout
      // edgeLogits    : pre-sigmoid logits emitted by the ONNX model (size E)
      // hitEnergiesMeV: per-node raw energies in MeV (size N)
      //
      // Returns a vector of length N: labels[i] = cluster ID (>= 0) or
      // -1 (unclustered after min_hits / min_energy_mev cleanup).
      std::vector<int> assemble(int nNodes,
                                const std::vector<int64_t>& edgeIndex,
                                const std::vector<float>& edgeLogits,
                                const std::vector<float>& hitEnergiesMeV) const;

    private:
      Config cfg_;
  };

}

#endif
