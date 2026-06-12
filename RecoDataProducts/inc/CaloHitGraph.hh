#ifndef RecoDataProducts_CaloHitGraph_hh
#define RecoDataProducts_CaloHitGraph_hh
//
// Per-disk graph carrying the inputs the GNN clustering ONNX model
// expects, plus per-node back-references to the source CaloHits.
//
// Transient data product: emitted by CaloHitGraphMaker, consumed by
// CaloClusterMakerGNN in the same job. Not registered in
// classes_def.xml (no ROOT serialisation).
//
// Tensor layout matches the ONNX model interface in
// calorimeter/GNN/docs/onnx_deployment.md (§2 / §7):
//
//   x          : nNodes * 6 floats, row-major (one row per hit)
//   edgeIndex  : 2 * nEdges int64s, row-major (src row first, dst row second)
//   edgeAttr   : nEdges * 8 floats, row-major (one row per directed edge)
//
// Feature column order is fixed:
//   nodes  : log_e, t, x, y, r, e_rel
//   edges  : dx, dy, d, dt, dlog_e, asym_e, logsum_e, dr
// The cluster module asserts these names against the loaded model's
// metadata_props (16j handshake).
//

#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "canvas/Persistency/Common/Ptr.h"

#include <cstdint>
#include <vector>

namespace mu2e {

  struct CaloHitGraph
  {
    int nNodes  = 0;
    int nEdges  = 0;
    int diskID  = -1;

    std::vector<float>   x;          // size 6 * nNodes
    std::vector<int64_t> edgeIndex;  // size 2 * nEdges
    std::vector<float>   edgeAttr;   // size 8 * nEdges

    std::vector<art::Ptr<CaloHit>> caloHitPtrs;  // size nNodes
  };

  using CaloHitGraphCollection = std::vector<CaloHitGraph>;

}

#endif
