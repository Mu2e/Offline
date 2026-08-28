// VDResamplerConfigure_module.cc
// Makes one pass of the data set and determines which particles need training for the resampler,
// then generates ONE self-contained training fcl per (source, particle). "dataSourceTag" tags the
// data source and, together with the plan's versionTag, runNumber and the VD id, is embedded in
// every generated file name so several campaigns/sources/detectors coexist in one directory:
//   train fcl : <TrainModule>_<ver>_VD<id>_<src>_pdg<pdgTok>.fcl
//   model     : nts.mu2e.STMVDResamplerModel_VD<id>_pdg<pdgTok>_<role>.<ver>.<run6>_<src8>.dat
//               (see VDResamplerNameHelper)
// A breakdown of the per-particle hit counts is written to
//   etc.mu2e.STMVDResamplerConfigure_VD<id>_hitSummary.<ver>.<run6>_<src8>.txt
// which the generator (VDResamplerGenerateMix) parses to recover ver/id/src/run and rebuild the
// model names via the SAME shared helper. Per-particle training config (momentum basis, curriculum,
// ...) is resolved from the training plan; common_training_config supplies the shared VD geometry,
// versionTag, runNumber, and the training-source selection (trainingFromROOTFile).
// Yongyi Wu, Mar. 2026

// stdlib includes
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <ctime>
#include <algorithm>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GeneralUtilities/inc/ParameterSetFromFile.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/STMMC/inc/VDResamplerPtotResampler.hh" // parseStage1Method (validation)
#include "Offline/STMMC/inc/VDResamplerNameHelper.hh"     // summaryFileName / modelFileName / rootDumpFileName

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/exception.h"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

typedef unsigned long VolumeId_type;

namespace mu2e {
  namespace {
    // Fallback only: two-stage training needs enough events to fit both stages, so when a plan
    // entry does NOT set SBDMuseTwoStageTraining we default to two-stage iff nHits > this. When the
    // plan sets SBDMuseTwoStageTraining explicitly, that value wins and this is not consulted.
    constexpr int kTwoStageTrainingHitThreshold = 10000;

    // Directory prefix for a generated file name. An empty directory means "no prefix", so the
    // file lands in the current working directory — on the grid, the job's local sandbox. A
    // non-empty directory gets exactly one trailing '/', whether or not it already had one.
    inline std::string dirPrefix(const std::string& dir) {
      if (dir.empty()) return dir;
      return (dir.back() == '/') ? dir : dir + "/";
    }

    // fhicl key for a pdgId: negatives can't start with '-', so use 'm<abs>'
    // (matching the model-file token convention), positives use the bare number.
    inline std::string pdgPlanKey(int pdgId) {
      return "pdg_" + std::string(pdgId < 0 ? "m" : "") + std::to_string(std::abs(pdgId));
    }

    // Resolve common_training_config.versionTag for THIS job's data source.
    //
    // Accepts either form, so existing plans keep working:
    //   * a bare string  -- one campaign tag shared by every source (the original behaviour);
    //   * a sequence     -- one tag PER SOURCE, indexed by the source's position in
    //                       VDResampler::dataSourceNames() (EleBeam, MuBeam, TargetStops1809,
    //                       Neutrals). Sources are re-trained on their own schedules, so tying
    //                       every source to one tag forced a version bump on all of them
    //                       whenever any single one was regenerated.
    //
    // The sequence is indexed by the ENUMERATED source order, not by the order the entries
    // happen to be written in, so the list must be as long as dataSourceNames() and an
    // unenumerated dataSourceTag cannot use the per-source form (it has no index to look up).
    std::string resolveVersionTag(const fhicl::ParameterSet& commonConfig,
                                  const std::string& dataSourceTag,
                                  const std::string& trainingPlanFile)
    {
      const std::vector<std::string>& sources = VDResampler::dataSourceNames();

      // The enumerated source list, for the error messages below.
      auto sourceList = [&sources]() {
        std::ostringstream os;
        for (size_t i = 0; i < sources.size(); ++i) os << (i ? ", " : "") << sources[i];
        return os.str();
      };

      if (commonConfig.is_key_to_atom("versionTag"))
        return commonConfig.get<std::string>("versionTag");

      if (!commonConfig.is_key_to_sequence("versionTag"))
        throw cet::exception("VDResamplerConfigure")
          << "common_training_config.versionTag in " << trainingPlanFile << " must be either a "
          << "string (one tag for every source) or a sequence of strings (one tag per source, "
          << "in dataSourceNames() order: " << sourceList() << ").";

      const std::vector<std::string> tags = commonConfig.get<std::vector<std::string>>("versionTag");
      if (tags.size() != sources.size())
        throw cet::exception("VDResamplerConfigure")
          << "common_training_config.versionTag in " << trainingPlanFile << " has " << tags.size()
          << " entries but there are " << sources.size() << " enumerated data sources; the "
          << "per-source form needs exactly one tag per source, in this order: " << sourceList()
          << ".";

      const int index = VDResampler::dataSourceIndex(dataSourceTag);
      if (index < 0)
        throw cet::exception("VDResamplerConfigure")
          << "dataSourceTag '" << dataSourceTag << "' is not one of the enumerated data sources ("
          << sourceList() << "), so it has no position in the per-source versionTag list in "
          << trainingPlanFile << ". Use a single versionTag string for an unenumerated source.";

      return tags[static_cast<size_t>(index)];
    }

    // Resolve a per-particle training entry from the plan ParameterSet. pdgId is the
    // OUTER axis (a particle usually trains similarly across sources); source is the
    // finer override. Fallback order:
    //   particle_training_config.<pdg_key>.<source>  ->  particle_training_config.<pdg_key>.default  ->
    //   particle_training_config.default.default
    // Throws if no default chain exists (the plan file is required and must supply at
    // least particle_training_config.default.default).
    fhicl::ParameterSet resolvePlanEntry(const fhicl::ParameterSet& plan,
                                         const std::string& source, int pdgId,
                                         const std::string& moduleName)
    {
      const fhicl::ParameterSet particle_training_config =
        plan.get<fhicl::ParameterSet>("particle_training_config", fhicl::ParameterSet());
      const std::string pdgKey = pdgPlanKey(pdgId);

      fhicl::ParameterSet pdgSet;
      if (particle_training_config.has_key(pdgKey))         pdgSet = particle_training_config.get<fhicl::ParameterSet>(pdgKey);
      else if (particle_training_config.has_key("default")) pdgSet = particle_training_config.get<fhicl::ParameterSet>("default");
      else
        throw cet::exception(moduleName)
          << "trainingPlanFile: 'particle_training_config' has neither \"" << pdgKey
          << "\" nor a 'default' entry.";

      if (pdgSet.has_key(source))            return pdgSet.get<fhicl::ParameterSet>(source);
      if (pdgSet.has_key("default"))         return pdgSet.get<fhicl::ParameterSet>("default");
      throw cet::exception(moduleName)
        << "trainingPlanFile: particle \"" << pdgKey << "\" has neither source \""
        << source << "\" nor a 'default' entry.";
    }

    // Per-particle stage-1 method string (validated against the resampler enum).
    std::string resolveStage1Method(const fhicl::ParameterSet& plan,
                                    const std::string& source, int pdgId,
                                    const std::string& moduleName)
    {
      const fhicl::ParameterSet entry = resolvePlanEntry(plan, source, pdgId, moduleName);
      const std::string method = entry.get<std::string>("stage1Method", "DIFFUSION");
      (void) VDResampler::parseStage1Method(method, moduleName); // validate; throws on typo
      return method;
    }

    // Plan-entry keys that are NOT VDResamplerTrain(FromRoot) parameters and so must not be
    // copied into the generated fcl. Only stage1Method is a generate-side / summary field;
    // every other key in a resolved entry (including SBDMuseTwoStageTraining and SBDMmomentumBasis)
    // is a training parameter that passes straight through.
    bool isPassThroughTrainingKey(const std::string& key) {
      return key != "stage1Method";
    }

    // Preferred emission order for the pass-through training keys. fhicl's ParameterSet does not
    // preserve the plan file's textual key order (get_names() is sorted, and @table:: expansion
    // scrambles it), so we impose a fixed, readable grouping here: architecture, then schedule,
    // then target/loss, then curriculum arrays, then planner. Any key not listed is emitted after
    // these, in fhicl's sorted order, so newly added plan keys still pass through (just ungrouped).
    const std::vector<std::string>& trainingKeyOrder() {
      static const std::vector<std::string> order = {
        // --- model architecture ---
        "SBDMuseTwoStageTraining", "SBDMmomentumBasis", "SBDMpositionBasis",
        "SBDMhidden", "SBDMlayers",
        "SBDMtimeEmbeddingDim", "SBDMinputEmbeddingDims", "SBDMconditionEmbeddingDims",
        // --- noise schedule ---
        "SBDMnoiseSchedule", "SBDMbetaMin", "SBDMbetaMax", "SBDMcosineOffset",
        "SBDMlogSigMin", "SBDMlogSigMax",
        // --- prediction target / loss ---
        "SBDMpredictionTarget", "SBDMlossWeightPower",
        // --- dim-weight controller / EMA ---
        "SBDMuseDimWeightController", "SBDMdimWeightControllerEMADecay",
        "SBDMuseEMANetwork", "SBDMemaNetworkDecay", "SBDMpromoteEMA",
        // --- curriculum (per-phase arrays) ---
        "SBDMtrainingCurriculumEpochs", "SBDMtrainingCurriculumGradientClip",
        "SBDMtrainingCurriculumLearningRate", "SBDMtrainingCurriculumBiasLowSigma",
        "SBDMtrainingCurriculumTLowBound",
        "SBDMtrainingCurriculumTFocusLow", "SBDMtrainingCurriculumTFocusHigh",
        "SBDMtrainingCurriculumTFocusFraction", "SBDMtrainingCurriculumBatchSize",
        "SBDMtrainingCurriculumSamplesDrawnPerEpoch",
        "SBDMtrainingCurriculumUseDimWeightController", "SBDMtrainingCurriculumUsePeakSampling",
        "SBDMtrainingCurriculumPromoteEMA", "SBDMtrainingCurriculumMinDelta",
        // --- planner ---
        "SBDMautoCurriculumPlanner", "SBDMplannerSmoothWindow", "SBDMplannerPatience",
        "SBDMplannerMinEpochsPerPhase", "SBDMplannerMaxEpochsPerPhase"
      };
      return order;
    }

    // Render a ParameterSet once and split it into a { key -> "key: value" text block } map,
    // reusing fhicl's own serializer (to_indented_string) so every value type — bool/int/double/
    // string/sequence/table — is formatted exactly as fhicl would write it. Multi-line values
    // (sequences/tables) keep their newlines; bracket-depth counting is used ONLY to detect where
    // a multi-line value ends (it is not an indentation control). The map lets the caller emit the
    // blocks in trainingKeyOrder() regardless of fhicl's internal sorted key order.
    std::map<std::string, std::string> renderKeyBlocks(const fhicl::ParameterSet& pSet) {
      std::map<std::string, std::string> blocks;
      std::istringstream iss(pSet.to_indented_string(0));
      std::string line;
      while (std::getline(iss, line)) {
        // A new top-level key begins at a line "key: ...". Find its name (up to the first ':').
        const std::size_t colon = line.find(':');
        if (colon == std::string::npos) continue;
        const std::string key = line.substr(0, colon);
        std::string block = line;
        int depth = 0;
        for (char c : line) { if (c == '[' || c == '{') ++depth; else if (c == ']' || c == '}') --depth; }
        // Absorb continuation lines until brackets balance (multi-line sequence/table value).
        while (depth > 0 && std::getline(iss, line)) {
          block += "\n" + line;
          for (char c : line) { if (c == '[' || c == '{') ++depth; else if (c == ']' || c == '}') --depth; }
        }
        blocks[key] = block;
      }
      return blocks;
    }
  }

  class VDResamplerConfigure : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs")};
        fhicl::Atom<art::InputTag> SimParticlemvTag{Name("SimParticlemvTag"), Comment("Tag identifying the SimParticlemv")};
        // VirtualDetectorID (and VDz0/VDr) are NOT module parameters: they come from the training
        // plan's common_training_config, the single source of truth used both to filter hits here
        // and to configure the generated training fcls.
        fhicl::Atom<std::string> VDResamplerDir{Name("VDResamplerDir"), Comment("Directory to store the generated summary (.txt) and model files. Defaults to empty, meaning no directory is prefixed and files land in the current working directory — on the grid, the job's local sandbox."), ""};
        fhicl::Atom<std::string> fclDir{Name("fclDir"), Comment("Directory to store the generated fhicl files"), ""};
        fhicl::Atom<std::string> dataSourceTag{Name("dataSourceTag"), Comment("A tag to distinguish different data sources, will be appended to the generated fcl and csv files (EleBeam, MuBeam, TargetStops1809, or Neutrals)")};
        fhicl::Atom<int> trainingThreshold{Name("trainingThreshold"), Comment("Minimum number of hits for a particle type to be included in the training"), 100};
        fhicl::Atom<std::string> trainingPlanFile{Name("trainingPlanFile"), Comment("Required fhicl guideline file with per-(source,pdg) training configuration (e.g. stage1Method). See resolveStage1Method for the lookup order.")};
        fhicl::Atom<bool> doROOTDump{Name("doROOTDump"), Comment("Whether to dump the VD hit info into a ROOT file for debugging and analysis"), true};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit VDResamplerConfigure(const Parameters& conf);
      void analyze(const art::Event& e);
      void endJob();
    private:
      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      art::ProductToken<SimParticleCollection> SimParticlemvToken;
      std::string VDResamplerDir;
      std::string fclDir;
      std::string dataSourceTag;
      int trainingThreshold;
      std::string trainingPlanFile;
      bool doROOTDump;
      // Training plan, loaded once in the constructor so analyze() can filter on the plan's VD id.
      // common_training_config is the single source of truth for VirtualDetectorID / VDz0 / VDr.
      fhicl::ParameterSet trainingPlan;
      fhicl::ParameterSet commonConfig;
      VolumeId_type VirtualDetectorID = 0;
      double VDz0 = 0.0;
      double VDr = 0.0;
      bool trainingFromROOT = true;   // common_training_config.trainingFromROOTFile (required)
      std::string versionTag;         // this source's tag, from common_training_config.versionTag
                                      // (required; a shared string or a per-source sequence —
                                      // see resolveVersionTag); embedded in file names
      int runNumber = 0;              // common_training_config.runNumber (required); 6-digit field in file names
      GlobalConstantsHandle<ParticleDataList> pdt;
      int pdgId = 0;
      double x = 0.0, y = 0.0, z = 0.0, px = 0.0, py = 0.0, pz = 0.0, mass = 0.0, E = 0.0, time = 0.0;
      VolumeId_type virtualdetectorId = 0;
      // Rejected-HIT counters (not events): analyze() loops over StepPointMCs, so one art event
      // can contribute several. Split by reason — a hit on a different VD says nothing about the
      // flux on ours, while a backward-going hit on our VD does.
      int nRejectedHitsWrongVD = 0;   // virtualdetectorId != VirtualDetectorID
      int nRejectedHitsBackward = 0;  // right VD, but pz <= 0 (upstream-going)
      TTree* ttree = nullptr;
      std::map<int, int> pdgIds; // <id, count>
  };

  VDResamplerConfigure::VDResamplerConfigure(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag())),
    VDResamplerDir(conf().VDResamplerDir()),
    fclDir(conf().fclDir()),
    dataSourceTag(conf().dataSourceTag()),
    trainingThreshold(conf().trainingThreshold()),
    trainingPlanFile(conf().trainingPlanFile()),
    doROOTDump(conf().doROOTDump()) {

    // Load the (required) training-plan guideline file once, up front. common_training_config is
    // the single source of truth for the VD geometry (id/z0/r), used both to filter hits in
    // analyze() and to configure the generated training fcls in endJob(). All three are required.
    if (trainingPlanFile.empty())
        throw cet::exception("VDResamplerConfigure")
            << "trainingPlanFile is required but empty.";
    trainingPlan = ParameterSetFromFile(trainingPlanFile).pSet();
    if (!trainingPlan.has_key("common_training_config"))
        throw cet::exception("VDResamplerConfigure")
            << "trainingPlanFile " << trainingPlanFile << " has no 'common_training_config' table.";
    commonConfig = trainingPlan.get<fhicl::ParameterSet>("common_training_config");
    // The VD geometry is required in common_training_config (single source of truth). Report each
    // missing key by name instead of letting fhicl throw a generic get<> error.
    auto requireCommon = [&](const char* key) {
      if (!commonConfig.has_key(key))
        throw cet::exception("VDResamplerConfigure")
          << "common_training_config in " << trainingPlanFile << " is missing required key '" << key << "'.";
    };
    requireCommon("versionTag");
    requireCommon("runNumber");
    requireCommon("trainingFromROOTFile");
    requireCommon("VirtualDetectorID");
    requireCommon("VDz0");
    requireCommon("VDr");
    // versionTag may be one shared string or one-per-source sequence; either way this job
    // ends up with the single tag for ITS dataSourceTag, so everything downstream (file
    // names, the summary the generator parses back) is unchanged.
    versionTag        = resolveVersionTag(commonConfig, dataSourceTag, trainingPlanFile);
    runNumber         = commonConfig.get<int>("runNumber");
    trainingFromROOT  = commonConfig.get<bool>("trainingFromROOTFile");
    VirtualDetectorID = static_cast<VolumeId_type>(commonConfig.get<int>("VirtualDetectorID"));
    VDz0 = commonConfig.get<double>("VDz0");
    VDr  = commonConfig.get<double>("VDr");

    if (doROOTDump) {
        art::ServiceHandle<art::TFileService> tfs;
        ttree = tfs->make<TTree>( "ttree", "Virtual Detectors Hit Summary");
        ttree->Branch("time", &time, "time/D"); // ns
        ttree->Branch("virtualdetectorId", &virtualdetectorId, "virtualdetectorId/l");
        ttree->Branch("pdgId", &pdgId, "pdgId/I");
        ttree->Branch("x", &x, "x/D"); // mm
        ttree->Branch("y", &y, "y/D"); // mm
        ttree->Branch("z", &z, "z/D"); // mm
        ttree->Branch("px", &px, "px/D"); // MeV
        ttree->Branch("py", &py, "py/D"); // MeV
        ttree->Branch("pz", &pz, "pz/D"); // MeV
        ttree->Branch("E", &E, "E/D"); // MeV
    }
  };

  void VDResamplerConfigure::analyze(const art::Event& event) {
    // Get the data products from the event
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    if (StepPointMCs.empty())
      return;
    auto const& SimParticles = event.getProduct(SimParticlemvToken);
    if (SimParticles.empty())
      return;

    // Loop over all VD hits
    for (const StepPointMC& step : StepPointMCs) {
      // Get the associated particle
      const SimParticle& particle = SimParticles.at(step.trackId());
      pdgId = particle.pdgId();

      virtualdetectorId = step.virtualDetectorId();
      pz = step.momentum().z();

      if (doROOTDump) {
          time = step.time();
          x = step.position().x();
          y = step.position().y();
          z = step.position().z();
          px = step.momentum().x();
          py = step.momentum().y();
          mass = pdt->particle(pdgId).mass();
          E = std::sqrt(step.momentum().mag2()+mass*mass)-mass; // Subtract the rest mass
          if (E < 0)
            throw cet::exception("LogicError", "Energy is negative");

          ttree->Fill();
      }

      // Filter hits on the virtual detector ID and pz, counting the two rejection reasons
      // separately. The wrong-VD case is checked first, so a hit is attributed to exactly one
      // reason and the two counters sum to the total number of rejected hits.
      if (virtualdetectorId != VirtualDetectorID)
      {
        // mf::LogWarning("VDResamplerConfigure") << "Rejected hit (wrong VD)\n"
        //                                        << "PDG ID = " << pdgId << ", VDID = " << virtualdetectorId << ", z = " << step.position().z() << ", pz = " << pz;
        nRejectedHitsWrongVD += 1;
        continue;
      }
      if (pz <= 0)
      {
        // mf::LogWarning("VDResamplerConfigure") << "Rejected hit (backward-going)\n"
        //                                        << "PDG ID = " << pdgId << ", VDID = " << virtualdetectorId << ", z = " << step.position().z() << ", pz = " << pz;
        nRejectedHitsBackward += 1;
        continue;
      }

      // Count the number of hits for each particle type for the summary
      if (pdgIds.find(pdgId) != pdgIds.end())
        pdgIds[pdgId] += 1;
      else
        pdgIds.emplace(std::make_pair(pdgId, 1));
    };
    return;
  };

  void VDResamplerConfigure::endJob() {
    mf::LogInfo log("Virtual Detector Resampler Training Configuration Summary");
    log << "========= Particle Summary =========\n";
    for (auto part : pdgIds)
      log << "PDGID " << part.first << ": " << part.second << "\n";
    log << "====================================\n";
    // Hit-level counts, matching the per-particle counts above (analyze() loops over StepPointMCs,
    // so one art event can contribute several). Accepted + rejected = every hit examined.
    int nAcceptedHits = 0;
    for (const auto& part : pdgIds) nAcceptedHits += part.second;
    const int nRejectedHits = nRejectedHitsWrongVD + nRejectedHitsBackward;
    log << "Accepted hits (VD" << VirtualDetectorID << ", pz>0): " << nAcceptedHits << "\n";
    log << "Rejected hits: " << nRejectedHits
        << " (wrong VD: " << nRejectedHitsWrongVD
        << ", backward-going on VD" << VirtualDetectorID << ": " << nRejectedHitsBackward << ")\n";

    // store a summary of the number of hits for each particle type in a comma-separated .txt file
    // for reference, all particles are included, but the particle types with hits below the training
    // threshold are marked with an asterisk in the TrainingMode column to indicate they are NOT trained.
    // Summary name embeds versionTag + VD id + dataSourceTag + runNumber (see VDResamplerNameHelper).
    // The generator parses all four back out of this file name, so it MUST be built via the shared helper.
    std::string summaryFile =
      dirPrefix(VDResamplerDir) + VDResampler::summaryFileName(versionTag, VirtualDetectorID, dataSourceTag, runNumber);
    std::string fclFilePath = fclDir.empty() ? VDResamplerDir : fclDir;
    std::ofstream sumOutFile(summaryFile);
    if (!sumOutFile)
        throw cet::exception("VDResamplerConfigure::endJob") << "Cannot open file " << summaryFile;

    // The ROOT dump this job writes is named by services.TFileService.fileName (set at launch,
    // outside this module). Log the recommended, version-consistent name so it matches the models.
    mf::LogInfo("VDResamplerConfigure::endJob")
        << "Recommended TFileService.fileName for this job's ROOT dump: "
        << VDResampler::rootDumpFileName(versionTag, VirtualDetectorID, dataSourceTag, runNumber)
        << " (set services.TFileService.fileName to keep it version-consistent).";

    // trainingPlan / commonConfig were loaded (and their required keys validated) in the
    // constructor. common_training_config is the single source of truth for the VD geometry
    // (VirtualDetectorID / VDz0 / VDr) and trainingFromROOT, which selects the training path for
    // EVERY generated job:
    //   true  -> VDResamplerTrainFromRoot module + EmptyEvent source, reading InputRootFile.
    //   false -> VDResamplerTrain module + the ART data source (StepPointMCs/SimParticles tags).
    // ROOT-source keys (only written when trainingFromROOT). InputRootFile may be @nil or absent,
    // in which case the generated fcl gets a placeholder + reminder instead of a path. fhicl
    // reports a @nil value as a present atom, so has_key/is_key_to_atom cannot distinguish it from
    // a real path — only the string conversion can, and it throws. Let the get<> attempt itself be
    // the test and treat a failure as "not set".
    const std::string treeName = commonConfig.get<std::string>("TreeName", "VDResamplerTrainingSetup/ttree");
    std::string inputRootFile;
    bool inputRootFileSet = false;
    if (commonConfig.has_key("InputRootFile")) {
      try {
        inputRootFile = commonConfig.get<std::string>("InputRootFile");
        inputRootFileSet = !inputRootFile.empty();
      } catch (const fhicl::exception&) {
        inputRootFileSet = false;  // @nil (or any non-string value): treat as unset
      }
    }
    // ART-source tags (only written when !trainingFromROOT).
    const std::string stepTag = commonConfig.get<std::string>("StepPointMCsTag", "compressDetStepMCsSTM116:");
    const std::string simTag  = commonConfig.get<std::string>("SimParticlemvTag", "compressDetStepMCsSTM116:");

    // Fixed 4-column schema (every row has the same field count for easy parsing):
    //   PDGID, HitCount, TrainingMode (1=all-at-once, 2=two-stage, *=below threshold, not trained),
    //   Stage1Method (DIFFUSION/INVERSE_CDF/... for two-stage; '-' otherwise).
    sumOutFile << "PDGID,HitCount,TrainingMode,Stage1Method\n";

    // Sanitize the data-source tag once for the module / path / fcl identifiers (the model and
    // summary file names sanitize internally via the shared helper). Same rule as the helper.
    const std::string sanitizedDataSourceTag = VDResampler::sanitizeTag(dataSourceTag);
    const std::string vdStr = std::to_string(VirtualDetectorID);

    for (const auto& part : pdgIds) {
      const int pdg   = part.first;
      const int nHits = part.second;
      sumOutFile << pdg << "," << nHits;

      if (nHits < trainingThreshold) {
        // Below threshold: not trained. Keep the 4-column schema (TrainingMode='*', Stage1Method='-').
        sumOutFile << ",*,-\n";
        mf::LogWarning("VDResamplerConfigure::endJob") << "Particle type with PDGID " << pdg
            << " has only " << nHits << " hits, which is below the training threshold of "
            << trainingThreshold << ". This particle type will NOT be included in the training.";
        continue;
      }

      // Resolve this particle's plan entry (pdg.<source> -> pdg.default -> default.default). The
      // entry already carries the @table::-expanded SBDM* keys plus SBDMuseTwoStageTraining /
      // SBDMmomentumBasis / stage1Method.
      const fhicl::ParameterSet entry =
        resolvePlanEntry(trainingPlan, dataSourceTag, pdg, "VDResamplerConfigure::endJob");
      // The plan's explicit SBDMuseTwoStageTraining wins; only when it is absent do we fall back to
      // the hit-count rule (two-stage iff there are enough events to fit both stages).
      const bool useTwoStageTraining =
        entry.get<bool>("SBDMuseTwoStageTraining", nHits > kTwoStageTrainingHitThreshold);
      // stage1Method (the V2 two-stage stage-1 pTotal source) is only meaningful for a two-stage
      // model; resolve/validate it only then. Single-stage rows write "-" in the summary column.
      const std::string stage1Method = useTwoStageTraining
        ? resolveStage1Method(trainingPlan, dataSourceTag, pdg, "VDResamplerConfigure::endJob")
        : std::string("-");

      // Summary columns: pdgId, hitCount, stage-mode (1=all-at-once, 2=two-stage), stage1Method.
      // VDResamplerGenerateMix reads the stage1Method column (fields[3]) only for two-stage particles.
      sumOutFile << "," << (useTwoStageTraining ? "2" : "1") << "," << stage1Method << "\n";

      // One self-contained fcl per (source, particle). Negative pdgId uses 'm<abs>' in file names.
      const std::string pdgIdstr = (pdg < 0) ? "m" + std::to_string(-pdg) : std::to_string(pdg);
      const std::string tag = "VD" + vdStr + sanitizedDataSourceTag + "pdg" + pdgIdstr;
      const std::string trainModule = trainingFromROOT ? "VDResamplerTrainFromRoot" : "VDResamplerTrain";
      const std::string moduleName  = trainModule + tag;
      const std::string pathName    = "trainPathVD" + vdStr + "pdg" + pdgIdstr;
      // fcl file name mirrors the training module actually used (…TrainFromRoot_… vs …Train_…).
      const std::string fclFile = dirPrefix(fclFilePath) + trainModule + "_" + VDResampler::sanitizeTag(versionTag)
                                + "_VD" + vdStr + "_" + sanitizedDataSourceTag + "_pdg" + pdgIdstr + ".fcl";

      std::ofstream fclOutFile(fclFile);
      if (!fclOutFile)
        throw cet::exception("VDResamplerConfigure::endJob") << "Cannot open file " << fclFile;

      // ---- header + source + services ----
      fclOutFile << "# This fcl is generated by VDResamplerConfigure_module.cc\n"
                 << "# Training configuration for VD" << vdStr << " " << dataSourceTag
                 << " pdg " << pdg << " (" << nHits << " hits).\n\n"
                 << "#include \"Offline/fcl/standardServices.fcl\"\n"
                 << "#include \"Production/JobConfig/pileup/STM/prolog.fcl\"\n\n";
      if (trainingFromROOT) {
        fclOutFile << "process_name: VDResamplerTrainParticle\n\n"
                   << "# No ART input: VDResamplerTrainFromRoot reads InputRootFile directly.\n"
                   << "source : {\n  module_type : EmptyEvent\n}\n\n";
      } else {
        fclOutFile << "process_name: VDResamplerTrain\n\n"
                   << "source : @local::STMPileup.VDResamplerSource\n\n";
      }
      fclOutFile << "services : @local::Services.Sim\n\n";

      // ---- analyzer block ----
      const std::string ind = "      ";
      fclOutFile << "physics: {\n  analyzers : {\n"
                 << "    " << moduleName << " : {\n"
                 << ind << "module_type : " << trainModule << "\n";

      // Source keys, per training path.
      if (trainingFromROOT) {
        if (inputRootFileSet) {
          fclOutFile << ind << "InputRootFile : \"" << inputRootFile << "\"\n";
        } else {
          // Plan left InputRootFile as @nil / unset: emit a placeholder + in-file reminder, and log one.
          fclOutFile << ind << "InputRootFile : @nil # TODO set the training ROOT file (left @nil in the plan)\n";
          mf::LogWarning("VDResamplerConfigure::endJob")
            << "InputRootFile is @nil / unset in the training plan; the generated fcl " << fclFile
            << " has a placeholder that MUST be set before running.";
        }
        fclOutFile << ind << "TreeName : \"" << treeName << "\"\n";
      } else {
        fclOutFile << ind << "StepPointMCsTag : \"" << stepTag << "\"\n"
                   << ind << "SimParticlemvTag : \"" << simTag << "\"\n";
      }

      // Model output files (.dat, full double precision), named via the shared helper so the
      // generator reconstructs identical paths. Under a resampler stage-1 method the 1-D pTotal
      // model is NOT trained (the generator draws pTotal directly), so omit stage1's file — the
      // train module trains stage1 only when that path is non-empty. Always train stage2.
      const std::string modelDir = dirPrefix(VDResamplerDir);
      if (useTwoStageTraining) {
        if (stage1Method == "DIFFUSION")
          fclOutFile << ind << "SBDMstage1ModelFile : \""
                     << modelDir << VDResampler::modelFileName("stage1", versionTag, VirtualDetectorID, dataSourceTag, pdg, runNumber) << "\"\n";
        fclOutFile << ind << "SBDMstage2ModelFile : \""
                   << modelDir << VDResampler::modelFileName("stage2", versionTag, VirtualDetectorID, dataSourceTag, pdg, runNumber) << "\"\n";
      } else {
        fclOutFile << ind << "SBDMallAtOnceModelFile : \""
                   << modelDir << VDResampler::modelFileName("allAtOnce", versionTag, VirtualDetectorID, dataSourceTag, pdg, runNumber) << "\"\n";
      }

      // VD geometry — always from the plan (single source of truth).
      fclOutFile << ind << "VirtualDetectorID : " << VirtualDetectorID << "\n"
                 << ind << "VDz0 : " << VDz0 << "\n"
                 << ind << "VDr : "  << VDr  << "\n"
                 << ind << "pdgID : " << pdg << "\n"
                 << ind << "SBDMtrainingSize : " << nHits << "\n";

      // Pass through the plan's training keys (SBDM*), in a fixed readable grouping. stage1Method is
      // excluded (generate-side). Keys not in trainingKeyOrder() are emitted afterwards in fhicl's
      // sorted order so newly added plan keys still flow through.
      const std::map<std::string, std::string> blocks = renderKeyBlocks(entry);
      std::vector<std::string> emitted;
      for (const std::string& key : trainingKeyOrder()) {
        auto it = blocks.find(key);
        if (it == blocks.end() || !isPassThroughTrainingKey(key)) continue;
        fclOutFile << ind << it->second << "\n";
        emitted.push_back(key);
      }
      for (const auto& kv : blocks) {
        if (!isPassThroughTrainingKey(kv.first)) continue;
        if (std::find(emitted.begin(), emitted.end(), kv.first) != emitted.end()) continue;
        fclOutFile << ind << kv.second << "\n";
      }

      fclOutFile << "    }\n  }\n";

      // ---- path + end_paths + services overrides ----
      fclOutFile << "  " << pathName << " : [" << moduleName << "]\n"
                 << "  end_paths: [" << pathName << "]\n"
                 << "}\n\n";
      // Remove the per-category log line limit so full training logs are kept.
      fclOutFile << "services.message.destinations.log.categories.default.limit: -1\n";
      // Seed from the wall clock so re-generated jobs don't all share one fixed seed.
      fclOutFile << "services.SeedService.baseSeed : " << (static_cast<long>(std::time(nullptr)) % 900000000 + 1) << "\n";
    }

    return;
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::VDResamplerConfigure)
