//
// Implementation of GnnGraphBuilder.
// See header for design notes / Python reference.
//

#include "Offline/CaloCluster/inc/GnnGraphBuilder.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"

#include "cetlib_except/exception.h"
#include "nlohmann/json.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

  // Canonical feature names — must agree with the norm sidecar emitted
  // by calorimeter/GNN/scripts/export_norm_stats.py and with the
  // model_version-stamped metadata_props of the .onnx artifacts.
  const std::vector<std::string> kNodeFeatures{
    "log_e", "t", "x", "y", "r", "e_rel"
  };
  const std::vector<std::string> kEdgeFeatures{
    "dx", "dy", "d", "dt", "dlog_e", "asym_e", "logsum_e", "dr"
  };

}

namespace mu2e {

  GnnGraphBuilder::Stats
  GnnGraphBuilder::loadStatsFromJson(const std::string& jsonPath)
  {
    std::ifstream in(jsonPath);
    if (!in.is_open()) {
      throw cet::exception("GnnGraphBuilder")
        << "cannot open norm sidecar: " << jsonPath;
    }
    nlohmann::json j;
    try {
      in >> j;
    } catch (const std::exception& e) {
      throw cet::exception("GnnGraphBuilder")
        << "JSON parse error in " << jsonPath << ": " << e.what();
    }

    auto require = [&](const char* k) {
      if (!j.contains(k)) {
        throw cet::exception("GnnGraphBuilder")
          << "norm sidecar missing key '" << k << "': " << jsonPath;
      }
    };
    require("schema_version");
    require("node_features");
    require("edge_features");
    require("node_mean");
    require("node_std");
    require("edge_mean");
    require("edge_std");

    const int schema = j["schema_version"].get<int>();
    if (schema != 1) {
      throw cet::exception("GnnGraphBuilder")
        << "unsupported norm sidecar schema_version " << schema
        << " (expected 1): " << jsonPath;
    }

    auto nodeNames = j["node_features"].get<std::vector<std::string>>();
    auto edgeNames = j["edge_features"].get<std::vector<std::string>>();
    if (nodeNames != kNodeFeatures) {
      throw cet::exception("GnnGraphBuilder")
        << "node_features mismatch in " << jsonPath;
    }
    if (edgeNames != kEdgeFeatures) {
      throw cet::exception("GnnGraphBuilder")
        << "edge_features mismatch in " << jsonPath;
    }

    Stats s;
    s.nodeMean = j["node_mean"].get<std::vector<float>>();
    s.nodeStd  = j["node_std"].get<std::vector<float>>();
    s.edgeMean = j["edge_mean"].get<std::vector<float>>();
    s.edgeStd  = j["edge_std"].get<std::vector<float>>();

    if (s.nodeMean.size() != kNodeFeatures.size()
        || s.nodeStd.size() != kNodeFeatures.size()) {
      throw cet::exception("GnnGraphBuilder")
        << "norm sidecar node stats size mismatch in " << jsonPath;
    }
    if (s.edgeMean.size() != kEdgeFeatures.size()
        || s.edgeStd.size() != kEdgeFeatures.size()) {
      throw cet::exception("GnnGraphBuilder")
        << "norm sidecar edge stats size mismatch in " << jsonPath;
    }
    for (float v : s.nodeStd) {
      if (v <= 0.0f) {
        throw cet::exception("GnnGraphBuilder")
          << "non-positive node_std entry in " << jsonPath;
      }
    }
    for (float v : s.edgeStd) {
      if (v <= 0.0f) {
        throw cet::exception("GnnGraphBuilder")
          << "non-positive edge_std entry in " << jsonPath;
      }
    }
    return s;
  }


  void GnnGraphBuilder::buildGraph(
      int diskID,
      const std::vector<const CaloHit*>& hits,
      const std::vector<art::Ptr<CaloHit>>& ptrs,
      const Calorimeter& cal,
      CaloHitGraph& out) const
  {
    const std::size_t n = hits.size();
    out.diskID  = diskID;
    out.nNodes  = static_cast<int>(n);
    out.caloHitPtrs = ptrs;
    out.x.assign(6 * n, 0.0f);
    out.edgeIndex.clear();
    out.edgeAttr.clear();
    out.nEdges = 0;
    if (n == 0) return;

    // Per-hit cached arrays: position, time, energy, radial.
    std::vector<double> xs(n), ys(n), ts(n), es(n), rs(n);
    double eMax = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const auto& h = *hits[i];
      const auto& cry = cal.crystal(h.crystalID());
      const auto& pos = cry.localPosition();
      xs[i] = pos.x();
      ys[i] = pos.y();
      ts[i] = h.time();
      es[i] = h.energyDep();
      rs[i] = std::sqrt(xs[i] * xs[i] + ys[i] * ys[i]);
      if (es[i] > eMax) eMax = es[i];
    }

    // Node features (6) + z-score normalisation.
    for (std::size_t i = 0; i < n; ++i) {
      const float feats[6] = {
        static_cast<float>(std::log1p(es[i])),
        static_cast<float>(ts[i]),
        static_cast<float>(xs[i]),
        static_cast<float>(ys[i]),
        static_cast<float>(rs[i]),
        static_cast<float>(eMax > 0.0 ? es[i] / eMax : 0.0)
      };
      for (int k = 0; k < 6; ++k) {
        out.x[6 * i + k] = (feats[k] - stats_.nodeMean[k]) / stats_.nodeStd[k];
      }
    }
    if (n == 1) return;

    // ----- 1. radius graph: brute-force pairwise i<j with r and dt cuts.
    const double rMax2 = cfg_.rMax * cfg_.rMax;
    std::vector<std::pair<int, int>> pairs;
    pairs.reserve(n * 4);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = i + 1; j < n; ++j) {
        const double dx = xs[i] - xs[j];
        const double dy = ys[i] - ys[j];
        if (dx * dx + dy * dy > rMax2) continue;
        if (std::abs(ts[i] - ts[j]) > cfg_.dtMax) continue;
        pairs.emplace_back(static_cast<int>(i), static_cast<int>(j));
      }
    }

    // Encode (s, d) pair as a single key for dedup. n is the count of
    // nodes; s and d are < n, so s*n + d uniquely identifies a directed
    // edge.
    auto encode = [n](int s, int d) -> long long {
      return static_cast<long long>(s) * static_cast<long long>(n)
           + static_cast<long long>(d);
    };

    std::unordered_set<long long> seen;
    std::vector<int> srcL, dstL;
    seen.reserve(4 * pairs.size() + 16);
    srcL.reserve(2 * pairs.size());
    dstL.reserve(2 * pairs.size());

    auto add_edge = [&](int s, int d) {
      const long long k = encode(s, d);
      if (seen.insert(k).second) {
        srcL.push_back(s);
        dstL.push_back(d);
      }
    };

    // Match the Python order: forward edges first (i->j for all pairs),
    // then backward edges (j->i). Keeps the pre-sort layout faithful
    // to the Python reference.
    for (const auto& p : pairs) add_edge(p.first, p.second);
    for (const auto& p : pairs) add_edge(p.second, p.first);

    // ----- 2. degree before kNN fallback.
    std::vector<int> degree(n, 0);
    for (int s : srcL) degree[s]++;

    // ----- 3. kNN fallback for under-connected nodes.
    for (std::size_t i = 0; i < n; ++i) {
      if (degree[i] >= static_cast<int>(cfg_.kMin)) continue;

      // Sorted candidate list (excluding self) by spatial distance.
      std::vector<std::pair<double, int>> cand;
      cand.reserve(n - 1);
      for (std::size_t j = 0; j < n; ++j) {
        if (j == i) continue;
        const double dx = xs[i] - xs[j];
        const double dy = ys[i] - ys[j];
        cand.emplace_back(std::sqrt(dx * dx + dy * dy),
                          static_cast<int>(j));
      }
      std::sort(cand.begin(), cand.end());

      // Mirror Python: tree.query returns k_query nearest including
      // self at j_pos=0. Self is silently skipped, so effectively
      // (k_query - 1) actual neighbours are inspected.
      const int kQuery = std::min(static_cast<int>(cfg_.kMin) * 3,
                                  static_cast<int>(n));
      const std::size_t maxN =
        std::min<std::size_t>(std::max(0, kQuery - 1), cand.size());
      int added = 0;
      for (std::size_t idx = 0; idx < maxN; ++idx) {
        const int j = cand[idx].second;
        if (std::abs(ts[i] - ts[static_cast<std::size_t>(j)]) > cfg_.dtMax)
          continue;
        add_edge(static_cast<int>(i), j);
        add_edge(j, static_cast<int>(i));
        ++added;
        if (degree[i] + added >= static_cast<int>(cfg_.kMin)) break;
      }
    }

    // ----- 4. Sort edges lexicographically by (src, dst) to match
    //         the Python deduplicate-by-encoded-value ordering.
    {
      std::vector<std::pair<int, int>> edges;
      edges.reserve(srcL.size());
      for (std::size_t k = 0; k < srcL.size(); ++k) {
        edges.emplace_back(srcL[k], dstL[k]);
      }
      std::sort(edges.begin(), edges.end());
      for (std::size_t k = 0; k < edges.size(); ++k) {
        srcL[k] = edges[k].first;
        dstL[k] = edges[k].second;
      }
    }

    // ----- 5. Per-source-node degree cap at k_max (keep nearest dsts).
    if (cfg_.kMax > 0) {
      std::vector<std::vector<int>> bySrc(n);
      for (std::size_t k = 0; k < srcL.size(); ++k) {
        bySrc[srcL[k]].push_back(static_cast<int>(k));
      }
      std::vector<bool> keep(srcL.size(), true);
      for (std::size_t s = 0; s < n; ++s) {
        if (bySrc[s].size() <= cfg_.kMax) continue;
        auto& idx = bySrc[s];
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
          const int da = dstL[a], db = dstL[b];
          const double dxa = xs[s] - xs[da];
          const double dya = ys[s] - ys[da];
          const double dxb = xs[s] - xs[db];
          const double dyb = ys[s] - ys[db];
          return (dxa * dxa + dya * dya) < (dxb * dxb + dyb * dyb);
        });
        for (std::size_t k = cfg_.kMax; k < idx.size(); ++k) {
          keep[idx[k]] = false;
        }
      }
      std::vector<int> sl, dl;
      sl.reserve(srcL.size());
      dl.reserve(srcL.size());
      for (std::size_t k = 0; k < srcL.size(); ++k) {
        if (keep[k]) {
          sl.push_back(srcL[k]);
          dl.push_back(dstL[k]);
        }
      }
      srcL.swap(sl);
      dstL.swap(dl);
    }

    // ----- 6. Build flat edge_index (src row first, dst row second).
    const int E = static_cast<int>(srcL.size());
    out.nEdges = E;
    out.edgeIndex.resize(2 * E);
    for (int e = 0; e < E; ++e) {
      out.edgeIndex[e]     = srcL[e];
      out.edgeIndex[E + e] = dstL[e];
    }

    // ----- 7. Edge features (8) + z-score normalisation.
    out.edgeAttr.resize(8 * E);
    for (int e = 0; e < E; ++e) {
      const int s = srcL[e];
      const int d = dstL[e];
      const double dx = xs[s] - xs[d];
      const double dy = ys[s] - ys[d];
      const double dist = std::sqrt(dx * dx + dy * dy);
      const double dt = ts[s] - ts[d];
      const double log_e_s = std::log1p(es[s]);
      const double log_e_d = std::log1p(es[d]);
      const double dlog_e = log_e_s - log_e_d;
      const double e_sum = es[s] + es[d];
      const double e_asym = (e_sum > 0.0)
                            ? (es[s] - es[d]) / e_sum : 0.0;
      const double logsum_e = std::log1p(e_sum);
      const double dr = rs[s] - rs[d];
      const double feats[8] = {dx, dy, dist, dt,
                               dlog_e, e_asym, logsum_e, dr};
      for (int k = 0; k < 8; ++k) {
        out.edgeAttr[8 * e + k] = static_cast<float>(
          (feats[k] - stats_.edgeMean[k]) / stats_.edgeStd[k]);
      }
    }
  }

}
