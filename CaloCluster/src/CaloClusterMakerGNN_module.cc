//
// CaloClusterMakerGNN — second half of the GNN clustering split design
// (calorimeter/GNN/docs/offline_integration.md §1.2).
//
// Constructor: load the ONNX session for one model artifact and assert
// the model's metadata_props match what the FHiCL config expects. Bail
// loudly on mismatch — silent tensor-layout drift after a retraining
// must not be possible.
//
// produce(): consumes CaloHitGraphCollection (from CaloHitGraphMaker),
// runs ONNX inference per disk, then assembles CaloClusters via the
// CCN+BFS10 recipe. The assembly logic (16f) is not yet implemented;
// produce() currently emits an empty CaloClusterCollection so the
// module compiles and links.
//
// The module is model-agnostic — the same C++ class instances run
// SimpleEdgeNet or CaloClusterNet (or any future model with the same
// tensor I/O signature). FHiCL parameters distinguish them:
//   - model_path
//   - expected_model_version
//   - expected_node_features / expected_edge_features
//   - tau_edge / bfs_expand_cut / min_hits / min_energy_mev
//   - output_instance
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/CaloCluster/inc/ClusterUtils.hh"
#include "Offline/CaloCluster/inc/GnnClusterAssembler.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloHitGraph.hh"

// Member-init order matters: env must outlive session_options, and
// session_options must outlive session. RAII resources lifetime
// captured in declaration order at the bottom of this class.
#include "onnxruntime/core/session/onnxruntime_cxx_api.h"

#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>


namespace {

  // Helper: split a comma-separated string into trimmed tokens.
  std::vector<std::string> splitCsv(const std::string& s)
  {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      // strip surrounding whitespace
      const auto a = tok.find_first_not_of(" \t");
      const auto b = tok.find_last_not_of(" \t");
      if (a == std::string::npos) continue;
      out.push_back(tok.substr(a, b - a + 1));
    }
    return out;
  }

  // Read an entry from the loaded model's metadata_props map by key.
  // Throws cet::exception if the key is absent.
  std::string readMetadataProp(const Ort::ModelMetadata& meta,
                               Ort::AllocatorWithDefaultOptions& alloc,
                               const char* key)
  {
    auto raw = meta.LookupCustomMetadataMapAllocated(key, alloc);
    if (!raw) {
      throw cet::exception("CaloClusterMakerGNN")
        << "ONNX metadata_props missing required key '" << key << "'";
    }
    return std::string(raw.get());
  }

}


namespace mu2e {

  class CaloClusterMakerGNN : public art::EDProducer
  {
    public:
      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> caloHitGraphCollection {
          Name("caloHitGraphCollection"),
          Comment("CaloHitGraphCollection input tag (from CaloHitGraphMaker)") };
        fhicl::Atom<std::string>   modelPath {
          Name("modelPath"),
          Comment("Relative path to the .onnx artifact; resolved by ConfigFileLookupPolicy") };
        fhicl::Atom<std::string>   expectedModelVersion {
          Name("expectedModelVersion"),
          Comment("Required value of metadata_props['model_version'] in the .onnx") };
        fhicl::Sequence<std::string> expectedNodeFeatures {
          Name("expectedNodeFeatures"),
          Comment("Required canonical node feature names; asserted against metadata_props") };
        fhicl::Sequence<std::string> expectedEdgeFeatures {
          Name("expectedEdgeFeatures"),
          Comment("Required canonical edge feature names; asserted against metadata_props") };
        fhicl::Atom<double>        tauEdge {
          Name("tauEdge"),
          Comment("Edge probability threshold (model-specific: 0.20 for CCN, 0.26 for SEN)") };
        fhicl::Atom<double>        bfsExpandCut {
          Name("bfsExpandCut"),
          Comment("BFS-style ExpandCut: hits below this energy join clusters but cannot recruit (MeV)"),
          10.0 };
        fhicl::Atom<unsigned>      minHits {
          Name("minHits"),
          Comment("Drop clusters with fewer than this many hits"),
          2 };
        fhicl::Atom<double>        minEnergyMeV {
          Name("minEnergyMeV"),
          Comment("Drop clusters with less than this much total energy (MeV)"),
          10.0 };
        fhicl::Atom<std::string>   outputInstance {
          Name("outputInstance"),
          Comment("Instance name on the emitted CaloClusterCollection (\"GNN\" by default)"),
          std::string("GNN") };
      };

      explicit CaloClusterMakerGNN(const art::EDProducer::Table<Config>& config);

      void produce(art::Event& event) override;

    private:
      art::ProductToken<CaloHitGraphCollection> graphToken_;

      std::string outputInstance_;
      std::vector<std::string> expectedNodeFeatures_;
      std::vector<std::string> expectedEdgeFeatures_;

      // CCN+BFS10 (or whatever the swappable model needs) recipe.
      GnnClusterAssembler assembler_;

      // ONNX Runtime resources. Member declaration order is intentional:
      // env_ must outlive sessionOptions_ which must outlive session_.
      Ort::Env                          env_;
      Ort::SessionOptions               sessionOptions_;
      Ort::Session                      session_;
      Ort::AllocatorWithDefaultOptions  allocator_;

      // Cached input/output names (returned by the session, owned by
      // RAII smart pointers from ONNX Runtime).
      std::vector<Ort::AllocatedStringPtr> inputNameHolders_;
      std::vector<Ort::AllocatedStringPtr> outputNameHolders_;
      std::vector<const char*>             inputNames_;
      std::vector<const char*>             outputNames_;
  };


  CaloClusterMakerGNN::CaloClusterMakerGNN(
      const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    graphToken_{consumes<CaloHitGraphCollection>(config().caloHitGraphCollection())},
    outputInstance_ (config().outputInstance()),
    expectedNodeFeatures_(config().expectedNodeFeatures()),
    expectedEdgeFeatures_(config().expectedEdgeFeatures()),
    assembler_({config().tauEdge(),
                config().bfsExpandCut(),
                config().minHits(),
                config().minEnergyMeV()}),
    env_(ORT_LOGGING_LEVEL_WARNING, "CaloClusterMakerGNN"),
    sessionOptions_(),
    session_{nullptr}
  {
    ConfigFileLookupPolicy lookup;
    const std::string modelPath = lookup(config().modelPath());

    // Move-construct the session in place (Ort::Session has no
    // copy-assign; the trick is to assign via std::move from a
    // temporary built with the resolved path).
    session_ = Ort::Session(env_, modelPath.c_str(), sessionOptions_);

    // Cache the input/output name strings — Run() takes const char*
    // arrays and we want them to live as long as the session does.
    {
      const std::size_t nIn  = session_.GetInputCount();
      const std::size_t nOut = session_.GetOutputCount();
      inputNameHolders_.reserve(nIn);
      outputNameHolders_.reserve(nOut);
      inputNames_.reserve(nIn);
      outputNames_.reserve(nOut);
      for (std::size_t i = 0; i < nIn; ++i) {
        inputNameHolders_.emplace_back(session_.GetInputNameAllocated(i, allocator_));
        inputNames_.push_back(inputNameHolders_.back().get());
      }
      for (std::size_t i = 0; i < nOut; ++i) {
        outputNameHolders_.emplace_back(session_.GetOutputNameAllocated(i, allocator_));
        outputNames_.push_back(outputNameHolders_.back().get());
      }
    }

    // Validate metadata_props against FHiCL expectations.
    Ort::ModelMetadata meta = session_.GetModelMetadata();
    const std::string gotVersion = readMetadataProp(meta, allocator_, "model_version");
    const std::string gotNodes   = readMetadataProp(meta, allocator_, "node_features");
    const std::string gotEdges   = readMetadataProp(meta, allocator_, "edge_features");

    if (gotVersion != config().expectedModelVersion()) {
      throw cet::exception("CaloClusterMakerGNN")
        << "model_version mismatch: ONNX has '" << gotVersion
        << "', FHiCL expected '" << config().expectedModelVersion()
        << "' (modelPath=" << modelPath << ")";
    }
    const auto gotNodeNames = splitCsv(gotNodes);
    if (gotNodeNames != expectedNodeFeatures_) {
      throw cet::exception("CaloClusterMakerGNN")
        << "node_features mismatch in " << modelPath
        << ": ONNX has '" << gotNodes
        << "', FHiCL expected " << expectedNodeFeatures_.size() << " entries";
    }
    const auto gotEdgeNames = splitCsv(gotEdges);
    if (gotEdgeNames != expectedEdgeFeatures_) {
      throw cet::exception("CaloClusterMakerGNN")
        << "edge_features mismatch in " << modelPath
        << ": ONNX has '" << gotEdges
        << "', FHiCL expected " << expectedEdgeFeatures_.size() << " entries";
    }

    produces<CaloClusterCollection>(outputInstance_);
  }


  void CaloClusterMakerGNN::produce(art::Event& event)
  {
    auto graphHandle = event.getHandle<CaloHitGraphCollection>(graphToken_);
    auto out = std::make_unique<CaloClusterCollection>();

    if (!graphHandle.isValid() || graphHandle->empty()) {
      event.put(std::move(out), outputInstance_);
      return;
    }

    const Calorimeter& cal = *(GeomHandle<Calorimeter>());
    const auto memInfo = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator,
                                                    OrtMemTypeDefault);

    for (const auto& graph : *graphHandle) {
      if (graph.nNodes <= 0 || graph.nEdges <= 0) continue;

      const int N = graph.nNodes;
      const int E = graph.nEdges;

      // Wrap the graph tensors as Ort::Value views — no copy.
      const std::array<int64_t, 2> xShape  {N, 6};
      const std::array<int64_t, 2> eiShape {2, E};
      const std::array<int64_t, 2> eaShape {E, 8};
      std::array<Ort::Value, 3> inputs{
        Ort::Value::CreateTensor<float>(memInfo,
            const_cast<float*>(graph.x.data()), graph.x.size(),
            xShape.data(), xShape.size()),
        Ort::Value::CreateTensor<int64_t>(memInfo,
            const_cast<int64_t*>(graph.edgeIndex.data()), graph.edgeIndex.size(),
            eiShape.data(), eiShape.size()),
        Ort::Value::CreateTensor<float>(memInfo,
            const_cast<float*>(graph.edgeAttr.data()), graph.edgeAttr.size(),
            eaShape.data(), eaShape.size())
      };

      auto outputs = session_.Run(Ort::RunOptions{nullptr},
                                  inputNames_.data(), inputs.data(), inputs.size(),
                                  outputNames_.data(), outputNames_.size());

      // Read edge_logits into a local buffer (don't mutate the
      // runtime's output).
      const Ort::Value& logitsVal = outputs.front();
      const float* logitsPtr = logitsVal.GetTensorData<float>();
      const std::size_t logitsCount =
        logitsVal.GetTensorTypeAndShapeInfo().GetElementCount();
      std::vector<float> edgeLogits(logitsPtr, logitsPtr + logitsCount);

      // Per-node energies (raw MeV) for the seed selection / expand
      // cut / min-energy cleanup. Pulled from the original CaloHits.
      std::vector<float> energiesMeV(N);
      for (int i = 0; i < N; ++i) {
        energiesMeV[i] = graph.caloHitPtrs[i]->energyDep();
      }

      // Run the CCN+BFS10 assembly on the directed-edge logits.
      const std::vector<int> labels = assembler_.assemble(
        N, graph.edgeIndex, edgeLogits, energiesMeV);

      // Bucket nodes by cluster ID.
      std::map<int, std::vector<int>> clusters;
      for (int i = 0; i < N; ++i) {
        if (labels[i] >= 0) clusters[labels[i]].push_back(i);
      }

      for (const auto& [cid, nodeIdxs] : clusters) {
        // Energy aggregates + seed identification.
        float energy = 0.0f, energyErrSq = 0.0f;
        int seedIdx = nodeIdxs.front();
        float seedEnergy = energiesMeV[seedIdx];
        CaloHitPtrVector caloHits;
        caloHits.reserve(nodeIdxs.size());
        for (int idx : nodeIdxs) {
          const auto& h = *graph.caloHitPtrs[idx];
          energy      += h.energyDep();
          energyErrSq += h.energyDepErr() * h.energyDepErr();
          caloHits.push_back(graph.caloHitPtrs[idx]);
          if (h.energyDep() > seedEnergy) {
            seedEnergy = h.energyDep();
            seedIdx    = idx;
          }
        }
        const auto& seedHit = *graph.caloHitPtrs[seedIdx];
        const float time    = seedHit.time();
        const float timeErr = seedHit.timeErr();

        // 3D centroid via the existing ClusterUtils helper (linear
        // energy weighting, matching the BFS module's convention).
        ClusterUtils utils(cal, caloHits, ClusterUtils::Linear);
        const CLHEP::Hep3Vector cog = utils.cog3Vector();

        out->emplace_back(graph.diskID, time, timeErr,
                          energy, std::sqrt(energyErrSq),
                          cog, caloHits,
                          static_cast<unsigned>(nodeIdxs.size()),
                          /*isSplit=*/false);
      }
    }

    event.put(std::move(out), outputInstance_);
  }

}

DEFINE_ART_MODULE(mu2e::CaloClusterMakerGNN)
