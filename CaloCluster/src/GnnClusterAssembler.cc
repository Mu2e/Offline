//
// Implementation of GnnClusterAssembler.
// See header for design notes / Python reference.
//

#include "Offline/CaloCluster/inc/GnnClusterAssembler.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <map>
#include <numeric>
#include <utility>
#include <vector>


namespace mu2e {

  std::vector<int>
  GnnClusterAssembler::assemble(
      int nNodes,
      const std::vector<int64_t>& edgeIndex,
      const std::vector<float>& edgeLogits,
      const std::vector<float>& hitEnergiesMeV) const
  {
    std::vector<int> labels(nNodes, -1);
    if (nNodes == 0) return labels;

    const int E = static_cast<int>(edgeLogits.size());
    // edgeIndex is flat (2E): src row first (indices 0..E-1) then dst
    // row (indices E..2E-1). Sanity check.
    if (static_cast<int>(edgeIndex.size()) != 2 * E) {
      // Defensive: malformed input; nothing to cluster.
      return labels;
    }

    // -------- 1. Sigmoid + symmetrise.
    // For each unordered pair {i, j}, accumulate the directed
    // probabilities; the threshold check uses their mean.
    std::map<std::pair<int, int>, std::pair<double, int>> pairScores;
    for (int e = 0; e < E; ++e) {
      const int s = static_cast<int>(edgeIndex[e]);
      const int d = static_cast<int>(edgeIndex[E + e]);
      const std::pair<int, int> key{std::min(s, d), std::max(s, d)};
      const double p = 1.0 / (1.0 + std::exp(-static_cast<double>(edgeLogits[e])));
      auto it = pairScores.find(key);
      if (it == pairScores.end()) {
        pairScores.emplace(key, std::make_pair(p, 1));
      } else {
        it->second.first  += p;
        it->second.second += 1;
      }
    }

    // -------- 2. Threshold + adjacency list.
    std::vector<std::vector<int>> adjList(nNodes);
    for (const auto& [key, sumCount] : pairScores) {
      const double avg = sumCount.first / static_cast<double>(sumCount.second);
      if (avg < cfg_.tauEdge) continue;
      adjList[key.first].push_back(key.second);
      adjList[key.second].push_back(key.first);
    }

    // -------- 3. BFS traversal with bfsExpandCut.
    // Seed selection: process nodes in descending energy order, so
    // each new cluster is rooted at its highest-energy hit (matching
    // Offline ClusterFinder semantics and the Python reference).
    std::vector<int> seedOrder(nNodes);
    std::iota(seedOrder.begin(), seedOrder.end(), 0);
    std::sort(seedOrder.begin(), seedOrder.end(),
              [&](int a, int b) {
                return hitEnergiesMeV[a] > hitEnergiesMeV[b];
              });

    int clusterId = 0;
    for (int seed : seedOrder) {
      if (labels[seed] >= 0) continue;
      std::deque<int> queue;
      queue.push_back(seed);
      labels[seed] = clusterId;
      while (!queue.empty()) {
        const int node = queue.front();
        queue.pop_front();
        // Hits below bfsExpandCut join the cluster but cannot recruit
        // further neighbours. They become leaves in the traversal.
        if (hitEnergiesMeV[node] < cfg_.bfsExpandCut) continue;
        for (int neigh : adjList[node]) {
          if (labels[neigh] < 0) {
            labels[neigh] = clusterId;
            queue.push_back(neigh);
          }
        }
      }
      ++clusterId;
    }

    // -------- 4. Cleanup: minHits.
    {
      std::map<int, int> hitCount;
      for (int label : labels) {
        if (label >= 0) hitCount[label]++;
      }
      for (int& label : labels) {
        if (label >= 0 && hitCount[label] < static_cast<int>(cfg_.minHits)) {
          label = -1;
        }
      }
    }

    // -------- 5. Cleanup: minEnergyMeV.
    {
      std::map<int, double> totalE;
      for (int i = 0; i < nNodes; ++i) {
        if (labels[i] >= 0) totalE[labels[i]] += hitEnergiesMeV[i];
      }
      for (int& label : labels) {
        if (label >= 0 && totalE[label] < cfg_.minEnergyMeV) {
          label = -1;
        }
      }
    }

    // -------- 6. Relabel to contiguous IDs in first-appearance order.
    {
      std::map<int, int> remap;
      int nextId = 0;
      for (int label : labels) {
        if (label >= 0 && remap.find(label) == remap.end()) {
          remap[label] = nextId++;
        }
      }
      for (int& label : labels) {
        if (label >= 0) label = remap[label];
      }
    }

    return labels;
  }

}
