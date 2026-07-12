// VDResamplerConfigure_module.cc
// Makes one pass of the data set and determines which particles will need training for the
// resampler, and generates a fcl file for training. "dataSourceTag" specifies the tag of the
// data source, and will be appended to the generated fcl file and then csv files to
// distinguish between different data sources. The generated fcl files will be stored and
// named "VDResamplerTrain_[dataSourceTag].fcl".
// On the other hand, a breakdown of the number of hits for each particle type will be stored
// in "VDResamplerConfigure_[dataSourceTag]_HitSummary.csv" for reference.
// Later when generating new samples using the resampler, the program will look for the
// generated fcl files to determine which particle source to use, which particle types to
// generate and which model parameters to load for each particle type.
// Yongyi Wu, Mar. 2026

// stdlib includes
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

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

#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

typedef unsigned long VolumeId_type;

namespace mu2e {
  namespace {
    constexpr int kTwoStageTrainingHitThreshold = 100000;

    // fhicl key for a pdgId: negatives can't start with '-', so use 'm<abs>'
    // (matching the model-file token convention), positives use the bare number.
    inline std::string pdgPlanKey(int pdgId) {
      return "pdg_" + std::string(pdgId < 0 ? "m" : "") + std::to_string(std::abs(pdgId));
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
  }

  class VDResamplerConfigure : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs")};
        fhicl::Atom<art::InputTag> SimParticlemvTag{Name("SimParticlemvTag"), Comment("Tag identifying the SimParticlemv")};
        fhicl::Atom<int> VirtualDetectorID{Name("VirtualDetectorID"), Comment("ID of the virtual detector to train on"), 116};
        fhicl::Atom<std::string> VDResamplerDir{Name("VDResamplerDir"), Comment("Directory to store the generated csv files")};
        fhicl::Atom<std::string> fclDir{Name("fclDir"), Comment("Directory to store the generated fhicl files"), ""};
        fhicl::Atom<std::string> dataSourceTag{Name("dataSourceTag"), Comment("A tag to distinguish different data sources, will be appended to the generated fcl and csv files (MuBeam, EleBeam, TgtStp1809, or Neutrals)")};
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
      VolumeId_type VirtualDetectorID = 0;
      std::string VDResamplerDir;
      std::string fclDir;
      std::string dataSourceTag;
      int trainingThreshold;
      std::string trainingPlanFile;
      bool doROOTDump;
      GlobalConstantsHandle<ParticleDataList> pdt;
      int pdgId = 0;
      double x = 0.0, y = 0.0, z = 0.0, px = 0.0, py = 0.0, pz = 0.0, mass = 0.0, E = 0.0, time = 0.0;
      VolumeId_type virtualdetectorId = 0;
      int nThrownEvents = 0;
      TTree* ttree = nullptr;
      std::map<int, int> pdgIds; // <id, count>
  };

  VDResamplerConfigure::VDResamplerConfigure(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag())),
    VirtualDetectorID(conf().VirtualDetectorID()),
    VDResamplerDir(conf().VDResamplerDir()),
    fclDir(conf().fclDir()),
    dataSourceTag(conf().dataSourceTag()),
    trainingThreshold(conf().trainingThreshold()),
    trainingPlanFile(conf().trainingPlanFile()),
    doROOTDump(conf().doROOTDump()) {
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

      if (virtualdetectorId != VirtualDetectorID || pz <= 0)
      {
        // mf::LogWarning("VDResamplerConfigure") << "Thrown event\n"
        //                                        << "PDG ID = " << pdgId << ", VDID = " << virtualdetectorId << ", z = " << step.position().z() << ", pz = " << pz;
        nThrownEvents += 1;
        continue; // Filter hits based on the virtual detector ID and pz
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
    log << "Number of thrown events: " << nThrownEvents << "\n";

    // store a summary of the number of hits for each particle type in a csv file for reference,
    // all particles are included, but the particle types with hits below the training threshold are
    // marked with an asterisk to indicate that they will not be included in the training.
    std::string summaryFile = VDResamplerDir + "/VDResamplerConfigure_" + dataSourceTag + "_HitSummary.csv";
    std::string fclFilePath = fclDir.empty() ? VDResamplerDir : fclDir;
    std::string fclFile = fclFilePath + "/VDResamplerTrain_" + dataSourceTag + ".fcl";
    std::ofstream sumOutFile(summaryFile);
    std::ofstream fclOutFile(fclFile);
    if (!sumOutFile) {
        throw cet::exception("VDResamplerConfigure::endJob") << "Cannot open file " << summaryFile;
    }
    if (!fclOutFile) {
        throw cet::exception("VDResamplerConfigure::endJob") << "Cannot open file " << fclFile;
    }

    // Load the (required) training-plan guideline file. Per-particle training choices
    // (stage1Method, and later curriculum settings) are resolved from it by (pdgId, source).
    if (trainingPlanFile.empty())
        throw cet::exception("VDResamplerConfigure::endJob")
            << "trainingPlanFile is required but empty.";
    const fhicl::ParameterSet trainingPlan =
        ParameterSetFromFile(trainingPlanFile).pSet();

    std::string trainingPathNames = "";
    std::vector<std::pair<std::string, std::string>> trainingPaths; // pairs of (moduleName, pathName)

    sumOutFile << "PDGID,HitCount,TrainingMode\n";

    fclOutFile << "# This fcl is generated by VDResamplerConfigure_module.cc\n";
    fclOutFile << "# Training configuration for the VD resampler is generated for each particle type\n\n";
    fclOutFile << "#include \"Offline/fcl/standardServices.fcl\"\n";
    fclOutFile << "#include \"Production/JobConfig/pileup/STM/prolog.fcl\"\n\n";
    fclOutFile << "process_name: VDResamplerTrain\n\n";
    fclOutFile << "source : @local::STMPileup.VDResamplerSource\n\n";
    fclOutFile << "services : @local::Services.Sim\n\n";
    fclOutFile << "physics: {\n  analyzers : {\n";

    for (const auto& part : pdgIds) {
      sumOutFile << part.first << "," << part.second;
      bool useTwoStageTraining = false;
      if (part.second < kTwoStageTrainingHitThreshold) {
        // if data is limited, use two-stage training to improve performance and reduce the required training data size for each model.
        useTwoStageTraining = true;
      }
      if (part.second < trainingThreshold) {
        sumOutFile << ",*\n";
        mf::LogWarning("VDResamplerConfigure::endJob") << "Particle type with PDGID " << part.first
            << " has only " << part.second << " hits, which is below the training threshold of "
            << trainingThreshold << ". This particle type will NOT be included in the training.";
      } else {
        // Columns: pdgId, hitCount, stage-mode (1=all-at-once, 2=two-stage), stage1Method.
        // stage1Method (the V2 two-stage stage-1 pTotal source: DIFFUSION / INVERSE_CDF /
        // SPLINE_CDF / KDE) is resolved per (pdgId, dataSourceTag) from the training plan.
        // VDResamplerGenerateMix reads this column (fields[3]) per particle.
        const std::string stage1Method =
          resolveStage1Method(trainingPlan, dataSourceTag, part.first, "VDResamplerConfigure::endJob");
        sumOutFile << "," << (useTwoStageTraining ? "2" : "1") << "," << stage1Method << "\n";

        // Generate the fcl file for training for this particle type
        // if pdgID is negative, we will use "m" instead of "-" in the filename to avoid issues with file naming
        std::string pdgIdstr = (part.first < 0) ? "m" + std::to_string(-part.first) : std::to_string(part.first);
        std::string sanitizedDataSourceTag = dataSourceTag;
        for (char& c : sanitizedDataSourceTag) {
          const bool isAlphaNum = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9');
          if (!isAlphaNum && c != '_') {
            c = '_';
          }
        }
        std::string moduleName = "VDResamplerTrainVD"+ std::to_string(VirtualDetectorID) + sanitizedDataSourceTag + "pdg" + pdgIdstr;
        std::string pathName = "trainPathVD" + std::to_string(VirtualDetectorID) + "pdg" + pdgIdstr;
        trainingPaths.push_back({moduleName, pathName});
        fclOutFile << "    " << moduleName << " : {\n"
                    << "      module_type : VDResamplerTrain\n"
                    << "      StepPointMCsTag : \"compressDetStepMCsSTM116:\"\n"
                    << "      SimParticlemvTag : \"compressDetStepMCsSTM116:\"\n"
                    << "      SBDMuseTwoStageTraining : " << (useTwoStageTraining ? "true" : "false") << "\n";
        if (useTwoStageTraining) {
          // Under a resampler stage-1 method the 1-D pTotal model is NOT trained (the
          // generator draws pTotal directly), so omit SBDMstage1ModelFile — the train
          // module trains stage1 only when that path is non-empty. Always train stage2.
          if (stage1Method == "DIFFUSION") {
            fclOutFile << "      SBDMstage1ModelFile : \"" << VDResamplerDir << "/SBDMmodel_stage1_VD" << VirtualDetectorID << "_" << sanitizedDataSourceTag << "_pdg" << pdgIdstr << ".csv\"\n";
          }
          fclOutFile << "      SBDMstage2ModelFile : \"" << VDResamplerDir << "/SBDMmodel_stage2_VD" << VirtualDetectorID << "_" << sanitizedDataSourceTag << "_pdg" << pdgIdstr << ".csv\"\n";
        } else {
          fclOutFile << "      SBDMallAtOnceModelFile : \"" << VDResamplerDir << "/SBDMmodel_allAtOnce_VD" << VirtualDetectorID << "_" << sanitizedDataSourceTag << "_pdg" << pdgIdstr << ".csv\"\n";
        }
        fclOutFile << "      VirtualDetectorID : " << VirtualDetectorID << "\n"
                    // << "      VDz0 : " << std::to_string(<value>) << "\n" // include if not VD116
                    // << "      VDr : " << std::to_string(<value>) << "\n" // include if not VD116
                    << "      pdgID : " << part.first << "\n"
                    // << "      SBDMhidden : 128\n"
                    // << "      SBDMlayers : 4\n"
                    // << "      SBDMoptimizer : \"ADAM\"\n"
                    // << "      SBDMadamBeta1 : 0.9\n"
                    // << "      SBDMadamBeta2 : 0.999\n"
                    // << "      SBDMadamEps : 1e-8\n"
                    // << "      SBDMnoiseSchedule : \"COSINE\"\n"
                    // << "      SBDMbetaMin : 1e-4\n"
                    // << "      SBDMbetaMax : 0.02\n"
                    // << "      SBDMcosineOffset : 0.008\n"
                    // << "      SBDMbatchSize : 32\n"
                    // << "      SBDMgradientClip : 1.0\n"
                    // << "      SBDMlearningRate : 1e-3\n"
                    // << "      SBDMdiffusionSteps : 200\n"
                    << "      SBDMtrainingSize : " << part.second << "\n"
                    // << "      SBDMtrainingEpochs : 10\n"
                    << "    }\n";
      }
    }
    fclOutFile << "  }\n";

    // Create trigger paths for each training module
    for (size_t i = 0; i < trainingPaths.size(); ++i) {
      fclOutFile << "  " << trainingPaths[i].second << " : [" << trainingPaths[i].first << "]\n";
      if (i > 0) trainingPathNames += ", ";
      trainingPathNames += trainingPaths[i].second;
    }

    fclOutFile << "  end_paths: [" + trainingPathNames + "]\n";
    fclOutFile << "}\n\n";
    fclOutFile << "services.SeedService.baseSeed : 8\n"; // fixed random seed for reproducibility, can be changed if needed

    return;
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::VDResamplerConfigure)
