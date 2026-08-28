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
#include <string>
#include <vector>
#include <map>

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
#include "Offline/STMMC/inc/VDResamplerConfigureCommon.hh" // plan resolution + fcl emission (shared)

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/exception.h"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

typedef unsigned long VolumeId_type;

namespace mu2e {
  // Plan resolution, the two-stage / stage1Method decision, and fcl emission now live in
  // VDResamplerConfigureCommon.hh so a re-emit job can regenerate the training fcls from an
  // existing hit summary without re-running this module's event loop. Writing the summary stays
  // here: it needs the hit counts, which only this pass over the data can produce.
  using VDResampler::EmitContext;

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
      // Everything the shared fcl emitter needs, resolved once in the constructor.
      EmitContext emitCtx;
      // VD id is mirrored out of emitCtx because analyze() filters hits on it; the rest of the
      // geometry / versioning lives in emitCtx and is read from there by the shared emitter.
      VolumeId_type VirtualDetectorID = 0;
      std::string versionTag;         // this source's tag; also used for the summary file name
      int runNumber = 0;              // 6-digit sequencer field in the summary file name
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
    // The VD geometry is required in common_training_config (single source of truth).
    VDResampler::requireCommonKeys(commonConfig, trainingPlanFile, "VDResamplerConfigure");
    // Everything the fcl emitter needs, resolved once. versionTag may be one shared string or a
    // one-per-source sequence; either way this job ends up with the single tag for ITS
    // dataSourceTag, so everything downstream (file names, the summary the generator parses back)
    // is unchanged.
    emitCtx = VDResampler::makeEmitContext(commonConfig, dataSourceTag, trainingPlanFile,
                                           VDResamplerDir, fclDir, "VDResamplerConfigure");
    // analyze() filters hits on the plan's VD id, and endJob() names the summary file from the
    // versioning fields; everything else the emitter needs stays in emitCtx.
    versionTag        = emitCtx.versionTag;
    runNumber         = emitCtx.runNumber;
    VirtualDetectorID = static_cast<VolumeId_type>(emitCtx.virtualDetectorID);

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
      VDResampler::detail::dirPrefix(VDResamplerDir)
      + VDResampler::summaryFileName(versionTag, VirtualDetectorID, dataSourceTag, runNumber);
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
    // constructor, and emitCtx was resolved from them there — it carries the VD geometry, the
    // versioning fields, and trainingFromROOT, which selects the training path for EVERY
    // generated job:
    //   true  -> VDResamplerTrainFromRoot module + EmptyEvent source, reading InputRootFile.
    //   false -> VDResamplerTrain module + the ART data source (StepPointMCs/SimParticles tags).

    // Fixed 4-column schema (every row has the same field count for easy parsing):
    //   PDGID, HitCount, TrainingMode (1=all-at-once, 2=two-stage, *=below threshold, not trained),
    //   Stage1Method (DIFFUSION/INVERSE_CDF/... for two-stage; '-' otherwise).
    sumOutFile << "PDGID,HitCount,TrainingMode,Stage1Method\n";

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
        VDResampler::resolvePlanEntry(trainingPlan, dataSourceTag, pdg, "VDResamplerConfigure::endJob");
      const bool useTwoStageTraining = VDResampler::resolveUseTwoStage(entry, nHits);
      // stage1Method (the V2 two-stage stage-1 pTotal source) is only meaningful for a two-stage
      // model; resolve/validate it only then. Single-stage rows write "-" in the summary column.
      const std::string stage1Method = useTwoStageTraining
        ? VDResampler::resolveStage1Method(trainingPlan, dataSourceTag, pdg, "VDResamplerConfigure::endJob")
        : VDResampler::kNoStage1Method();

      // Summary columns: pdgId, hitCount, stage-mode (1=all-at-once, 2=two-stage), stage1Method.
      // VDResamplerGenerateMix reads the stage1Method column (fields[3]) only for two-stage particles.
      sumOutFile << "," << (useTwoStageTraining ? "2" : "1") << "," << stage1Method << "\n";

      // One self-contained fcl per (source, particle), via the shared emitter.
      VDResampler::emitOneTrainingFcl(entry, emitCtx, pdg, nHits, useTwoStageTraining,
                                      stage1Method, "VDResamplerConfigure::endJob");
    }

    return;
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::VDResamplerConfigure)
