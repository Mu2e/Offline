// VDResamplerReconfigure_module.cc
// Re-emits the VDResampler training fcls from an EXISTING hit summary, without re-running over the
// data. VDResamplerConfigure must make a full pass over a dataset to count hits per particle; that
// pass is expensive and its result — the per-particle hit counts, plus the training-mode and
// stage-1-method decisions — is already recorded in
//   etc.mu2e.STMVDResamplerConfigure_VD<id>_hitSummary.<ver>.<run6>_<src8>.txt
// So when only the TRAINING PLAN changes (curriculum, batch schedule, architecture, ...), the fcls
// can be regenerated straight from that summary.
//
// The summary is strictly an INPUT here: this module never rewrites it. Only the hit-counting job
// has the information to produce one, and the recorded TrainingMode / Stage1Method are taken as
// given rather than re-derived, so a plan edit cannot silently flip an already-trained particle's
// stage mode and orphan its model files.
//
// versionTag, VD id, dataSourceTag and runNumber are recovered from the summary FILE NAME (the same
// four components VDResamplerGenerateMix parses back), so they cannot drift from the models the
// summary describes. Everything else — plan resolution, key ordering, file naming, the emitted text
// — is shared with VDResamplerConfigure via VDResamplerConfigureCommon.hh.
//
// No ART input: run it with an EmptyEvent source; all the work happens in endJob().
// Yongyi Wu, Aug. 2026

// stdlib includes
#include <string>
#include <vector>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/GeneralUtilities/inc/ParameterSetFromFile.hh"
#include "Offline/STMMC/inc/VDResamplerConfigureCommon.hh"
#include "Offline/STMMC/inc/VDResamplerNameHelper.hh"

namespace mu2e {

  class VDResamplerReconfigure : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<std::string> hitSummaryFile{Name("hitSummaryFile"),
          Comment("Existing hit summary (.txt) written by VDResamplerConfigure. Read-only: this "
                  "module never rewrites it. Its NAME supplies versionTag, VirtualDetectorID, "
                  "dataSourceTag and runNumber, so it must keep the "
                  "etc.mu2e.STMVDResamplerConfigure_VD<id>_hitSummary.<ver>.<run6>_<src8>.txt "
                  "convention.")};
        fhicl::Atom<std::string> trainingPlanFile{Name("trainingPlanFile"),
          Comment("Required fhicl guideline file with per-(source,pdg) training configuration. The "
                  "point of this module is to re-emit fcls after editing it.")};
        fhicl::Atom<std::string> VDResamplerDir{Name("VDResamplerDir"),
          Comment("Directory the generated fcls should point at for model (.dat) files. Must match "
                  "where the training jobs will write, since these paths are baked into the fcls. "
                  "Defaults to empty, meaning no prefix — the training job's working directory."), ""};
        fhicl::Atom<std::string> fclDir{Name("fclDir"),
          Comment("Directory to write the generated fhicl files into. Defaults to empty, which "
                  "falls back to VDResamplerDir (and, if that is empty too, to this job's working "
                  "directory). Set it when the models belong in a shared area but the fcls should "
                  "land locally."), ""};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit VDResamplerReconfigure(const Parameters& conf);
      void analyze(const art::Event&) override {}   // no ART input; everything happens in endJob()
      void endJob() override;
    private:
      std::string hitSummaryFile;
      std::string trainingPlanFile;
      std::string VDResamplerDir;
      std::string fclDir;
      fhicl::ParameterSet trainingPlan;
      fhicl::ParameterSet commonConfig;
  };

  VDResamplerReconfigure::VDResamplerReconfigure(const Parameters& conf) :
    art::EDAnalyzer(conf),
    hitSummaryFile(conf().hitSummaryFile()),
    trainingPlanFile(conf().trainingPlanFile()),
    VDResamplerDir(conf().VDResamplerDir()),
    fclDir(conf().fclDir()) {

    // Load and validate the training plan up front, so a bad plan fails before any file is written.
    if (trainingPlanFile.empty())
      throw cet::exception("VDResamplerReconfigure")
        << "trainingPlanFile is required but empty.";
    trainingPlan = ParameterSetFromFile(trainingPlanFile).pSet();
    if (!trainingPlan.has_key("common_training_config"))
      throw cet::exception("VDResamplerReconfigure")
        << "trainingPlanFile " << trainingPlanFile << " has no 'common_training_config' table.";
    commonConfig = trainingPlan.get<fhicl::ParameterSet>("common_training_config");
    VDResampler::requireCommonKeys(commonConfig, trainingPlanFile, "VDResamplerReconfigure");
  }

  void VDResamplerReconfigure::endJob() {
    // Recover (versionTag, VD id, dataSourceTag, runNumber) from the summary FILE NAME — the same
    // four components VDResamplerGenerateMix parses back out of it. Taking them from the name
    // rather than from the plan is what keeps the regenerated fcls pointing at the very models this
    // summary describes, even if the plan's versionTag/runNumber have since moved on.
    std::string versionTag, dataSourceTag;
    int summaryVDID = 0;
    int runNumber = 0;
    if (!VDResampler::parseSummaryFileName(hitSummaryFile, versionTag, summaryVDID, dataSourceTag, runNumber))
      throw cet::exception("VDResamplerReconfigure")
        << "Hit summary file name does not match the expected "
        << "etc.mu2e.STMVDResamplerConfigure_VD<id>_hitSummary.<ver>.<run6>_<src8>.txt convention: "
        << hitSummaryFile;

    const std::vector<VDResampler::SummaryRow> rows =
      VDResampler::readHitSummary(hitSummaryFile, "VDResamplerReconfigure");

    // Build the emit context from the plan, then override the four naming components with the ones
    // recovered above. makeEmitContext() would otherwise take versionTag/runNumber from the plan,
    // which is exactly the drift we are avoiding.
    VDResampler::EmitContext ctx =
      VDResampler::makeEmitContext(commonConfig, dataSourceTag, trainingPlanFile,
                                   VDResamplerDir, fclDir, "VDResamplerReconfigure");
    ctx.versionTag  = versionTag;
    ctx.runNumber   = runNumber;
    ctx.generatedBy = "VDResamplerReconfigure_module.cc (from " + hitSummaryFile + ")";

    // The VD id in the summary name must match the plan's, otherwise the plan being applied was
    // written for a different detector than the models this summary points at.
    if (static_cast<unsigned long>(summaryVDID) != ctx.virtualDetectorID)
      throw cet::exception("VDResamplerReconfigure")
        << "Summary file " << hitSummaryFile << " is for VD" << summaryVDID
        << " but common_training_config.VirtualDetectorID in " << trainingPlanFile
        << " is " << ctx.virtualDetectorID << ".";

    mf::LogInfo("VDResamplerReconfigure")
        << "Re-emitting training fcls for VD" << summaryVDID << " " << dataSourceTag
        << " (version " << versionTag << ", run " << runNumber << ") from " << hitSummaryFile
        << ": " << rows.size() << " trained particle(s), plan " << trainingPlanFile << ".";

    for (const VDResampler::SummaryRow& row : rows) {
      // Resolve this particle's plan entry (pdg.<source> -> pdg.default -> default.default), then
      // emit with the stage mode / stage-1 method AS RECORDED in the summary.
      const fhicl::ParameterSet entry =
        VDResampler::resolvePlanEntry(trainingPlan, dataSourceTag, row.pdgId, "VDResamplerReconfigure");

      // A plan edit that changes a particle's stage mode cannot be honoured here: the summary and
      // the existing model files were built the other way. Warn rather than throw, since the rest
      // of the plan's changes are still applied correctly.
      if (entry.has_key("SBDMuseTwoStageTraining")
          && entry.get<bool>("SBDMuseTwoStageTraining") != row.useTwoStage)
        mf::LogWarning("VDResamplerReconfigure")
          << "Plan sets SBDMuseTwoStageTraining=" << entry.get<bool>("SBDMuseTwoStageTraining")
          << " for pdg " << row.pdgId << ", but " << hitSummaryFile << " records it as "
          << (row.useTwoStage ? "two-stage" : "single-stage")
          << ". Keeping the summary's value so the generated fcl stays consistent with the model "
          << "file names; re-run VDResamplerConfigure over the data to change the stage mode.";

      const std::string fclFile =
        VDResampler::emitOneTrainingFcl(entry, ctx, row.pdgId, row.hitCount,
                                        row.useTwoStage, row.stage1Method, "VDResamplerReconfigure");
      mf::LogInfo("VDResamplerReconfigure")
        << "  pdg " << row.pdgId << " (" << row.hitCount << " hits, "
        << (row.useTwoStage ? "two-stage" : "all-at-once") << "): " << fclFile;
    }
  }

}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::VDResamplerReconfigure)
