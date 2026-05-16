//
// testGnnClusterAssembler_main — Stage-2 parity test for the C++
// GnnClusterAssembler (Task 16g).
//
// Reads a JSON parity payload (produced by
// calorimeter/GNN/scripts/dump_parity_payloads.py), replays
// GnnClusterAssembler::assemble for each graph, and asserts the
// emitted cluster_labels are byte-identical to the Python reference
// labels stored in the same payload. Exits non-zero on any mismatch.
//
// Usage:
//   testGnnClusterAssembler  <path-to-parity.json>
//
// Default path (relative to MUSE_WORK_DIR):
//   ../projects/calorimeter/GNN/tests/parity/calo_cluster_net_v2_stage1.parity.json
//

#include "Offline/CaloCluster/inc/GnnClusterAssembler.hh"

#include "nlohmann/json.hpp"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>


namespace {

  constexpr const char* kDefaultPayload =
    "/exp/mu2e/app/users/wzhou2/projects/calorimeter/GNN/tests/parity/calo_cluster_net_v2_stage1.parity.json";

}


int main(int argc, char* argv[])
{
  const std::string path = (argc > 1) ? argv[1] : kDefaultPayload;

  std::ifstream in(path);
  if (!in.is_open()) {
    std::cerr << "[FAIL] cannot open " << path << "\n";
    return 2;
  }
  nlohmann::json j;
  try { in >> j; }
  catch (const std::exception& e) {
    std::cerr << "[FAIL] JSON parse error in " << path
              << ": " << e.what() << "\n";
    return 2;
  }

  if (j.value("schema_version", 0) != 1) {
    std::cerr << "[FAIL] unsupported schema_version "
              << j.value("schema_version", 0) << "\n";
    return 2;
  }

  mu2e::GnnClusterAssembler::Config cfg;
  cfg.tauEdge      = j.at("tau_edge").get<double>();
  cfg.bfsExpandCut = j.at("bfs_expand_cut").get<double>();
  cfg.minHits      = j.at("min_hits").get<unsigned>();
  cfg.minEnergyMeV = j.at("min_energy_mev").get<double>();
  mu2e::GnnClusterAssembler asm_(cfg);

  const auto& graphs = j.at("graphs");
  std::cout << "Loaded " << graphs.size() << " graphs from " << path << "\n";
  std::cout << "Recipe: tauEdge=" << cfg.tauEdge
            << " bfsExpandCut=" << cfg.bfsExpandCut
            << " minHits=" << cfg.minHits
            << " minEnergyMeV=" << cfg.minEnergyMeV << "\n";
  std::cout << "model_version="
            << j.value("model_version", std::string("?")) << "\n";

  std::size_t nGraphs = 0;
  std::size_t nMismatchGraphs = 0;
  std::size_t nMismatchNodes = 0;
  std::size_t maxMismatchPerGraph = 0;

  for (const auto& g : graphs) {
    const int N = g.at("n_nodes").get<int>();
    const int E = g.at("n_edges").get<int>();
    auto edgeIndex   = g.at("edge_index").get<std::vector<int64_t>>();
    auto edgeLogits  = g.at("edge_logits").get<std::vector<float>>();
    auto energies    = g.at("energies").get<std::vector<float>>();
    auto refLabels   = g.at("cluster_labels").get<std::vector<int>>();

    if (static_cast<int>(edgeLogits.size()) != E
        || static_cast<int>(energies.size()) != N
        || static_cast<int>(edgeIndex.size()) != 2 * E
        || static_cast<int>(refLabels.size()) != N) {
      std::cerr << "[FAIL] graph " << nGraphs
                << ": shape mismatch in payload\n";
      return 2;
    }

    auto outLabels = asm_.assemble(N, edgeIndex, edgeLogits, energies);
    if (static_cast<int>(outLabels.size()) != N) {
      std::cerr << "[FAIL] graph " << nGraphs
                << ": assemble returned size " << outLabels.size()
                << " (expected " << N << ")\n";
      return 2;
    }

    std::size_t diffs = 0;
    for (int i = 0; i < N; ++i) {
      if (outLabels[i] != refLabels[i]) ++diffs;
    }
    if (diffs > 0) {
      ++nMismatchGraphs;
      nMismatchNodes += diffs;
      if (diffs > maxMismatchPerGraph) maxMismatchPerGraph = diffs;
      if (nMismatchGraphs <= 5) {
        std::cerr << "[diff] graph " << nGraphs
                  << " (N=" << N << ", E=" << E
                  << "): " << diffs << " mismatched node(s)\n";
        std::cerr << "       py : ";
        for (int v : refLabels) std::cerr << v << " ";
        std::cerr << "\n       cpp: ";
        for (int v : outLabels) std::cerr << v << " ";
        std::cerr << "\n";
      }
    }
    ++nGraphs;
  }

  std::cout << "\n=== Summary ===\n";
  std::cout << "graphs:         " << nGraphs << "\n";
  std::cout << "mismatch graphs: " << nMismatchGraphs << "\n";
  std::cout << "mismatch nodes:  " << nMismatchNodes << "\n";
  std::cout << "max diffs/graph: " << maxMismatchPerGraph << "\n";

  if (nMismatchGraphs > 0) {
    std::cerr << "[FAIL] cluster-label parity broken on "
              << nMismatchGraphs << " / " << nGraphs << " graph(s)\n";
    return 1;
  }
  std::cout << "[PASS] all " << nGraphs
            << " graphs match Python cluster_labels byte-exactly\n";
  return 0;
}
