// VDResamplerRegeneratePlots_module.cc
//
// Rebuild the VDResampler validation plots from an ALREADY-GENERATED ROOT dump tree,
// instead of from a live generation job.
//
// Why this exists
// ---------------
// A grid campaign runs VDResamplerToCrystals.fcl (or VDResamplerGenerateMix.fcl) as N
// independent jobs, each writing its own TFileService file holding the generated-sample
// dump tree AND a "validation/" directory of comparison plots. The dump TREES concatenate
// correctly under hadd, but the validation DIRECTORY does not:
//   * every job re-read the SAME mother file, so hadd sums N copies of the mother;
//   * both sides were already normalized to unit area per job, so summing them averages
//     jobs with unequal statistics as if they were equal;
//   * each job chose its rebin factor from ITS OWN sample count (the target is
//     ~2*sqrt(min(nMother,nGenerated))), so the merged statistics are stuck at
//     single-job resolution;
//   * the stored 2D z-ranges and the c_* / c2_* canvases are frozen per job.
// So the correct merge is: hadd the trees, DISCARD the per-job validation, and recompute
// the whole comparison once from the merged tree. That is what this module does.
//
// It is equally useful outside the grid case: point it at any single generated ROOT file
// to redraw its plots — after a change to the plotting code, against a different mother
// file, or simply to get the canvases back once a file has been stripped down to its tree.
//
// ONE SOURCE PER JOB
// ------------------
// This module deliberately handles a SINGLE data source at a time: hadd the jobs belonging
// to one source, run this on the result, repeat per source. A source fixes the sourceTag,
// hence one summary file, one mother file, and one plan lookup axis — so the only thing
// left to separate is the PARTICLE.
//
// That restriction is what keeps the configuration honest. The transform a set is plotted
// in is a property of the models that produced it, and those differ per (source, particle);
// mixing sources into one job would mean resolving a different basis per summaryIndex and
// silently plotting several coordinate systems into one directory tree. Fixing the source
// removes that axis entirely: every set in a run of this module shares the source, the
// mother file and the geometry, and differs only in pdgId — which the plan resolves on its
// own axis anyway.
//
// The tree's summaryIndex branch is therefore not used to partition anything. It is checked
// instead: a tree carrying more than one distinct value was produced by a multi-source mix
// and is rejected, since there is no single mother file such a sample could be compared
// against. See the SourceIndex parameter for the escape hatch.
//
// What it does
// ------------
// endJob() reads the generated dump tree, splits its entries by pdgId, books one
// ValidationPlots per particle in its own "validation/<label>/" subdirectory, refills the
// generated side entry by entry, and calls the same finalize() the generator calls — so the
// rebinning, normalization, metrics, canvases AND the per-particle validationMetrics tree
// are all produced by exactly the code path that produced the per-job ones, just with the
// full merged statistics behind them. Because ValidationPlots books everything (the
// histograms, the canvases and the metric tree alike) into the subdirectory it is handed,
// one particle's objects can never land in another's.
//
// The plots are written to this module's OWN TFileService directory, so the module never
// modifies its input: "regenerate" means produce a fresh file holding only the recomputed
// comparison, not edit the old one in place. Point services.TFileService.fileName at a new
// file and keep what this writes.
//
// The training plan is the single source of truth
// -----------------------------------------------
// Nothing about the transform or the geometry is restated in this module's fcl:
//   * hitSummaryFile's FILE NAME carries sourceTag/versionTag/runNumber, recovered with the
//     shared parseSummaryFileName — the same recovery loadSourceSummary does in the
//     generator.
//   * trainingPlanFile is then consulted per (sourceTag, pdgId) through the shared
//     resolvePlanEntry, whose (pdg -> source -> default) fallback is the one the plan file
//     documents. That is where each particle's SBDMmomentumBasis / SBDMpositionBasis comes
//     from, so a species trained in V3_..._TIME_ASINH is plotted in it while its neighbour
//     in V2_PTOT_SLOPES_ASINH is plotted in that.
//   * common_training_config supplies VirtualDetectorID / VDz0 / VDr, so the geometry
//     cannot drift from what training used.
// The mother distribution still needs the source ROOT file, which the summary does not
// name; that is the one thing listed here (MotherRootFile).
//
// The transformed-coordinate caveat
// ---------------------------------
// A live generation job hands fillGenerated() the model's OWN transformed output `g`
// alongside its inversion, which is what lets the transformed histograms show the model's
// output rather than a re-forward-transform of the inverse. The dump tree stores only the
// PHYSICAL values, so replaying from it must forward-transform them back. The transformed
// plots produced here are therefore a round trip (inverse then forward) of what the live
// job plotted directly.
//
// In practice this is close to a no-op — the transform pair is analytic and mutually
// inverse — but it is not bit-identical, and it CANNOT show an inversion asymmetry the way
// the live plots can, because any such asymmetry cancels in the round trip. The mother side
// goes through the very same forward transform (see ValidationPlots::fillMother), so the
// two sides are still treated identically and the comparison stays fair; it is only the
// live-job property of catching inversion bugs that is lost. Read the PHYSICAL plots as
// authoritative here.
//
// rawPtot is deliberately left invalid on the replayed GeneratedTransformed: it is a
// generation-path detail (the V2 resampler's raw draw) that the dump does not record, and
// ValidationPlots does not read it.
//
// Example fcl
// -----------
//   physics.analyzers.regenerate : {
//     module_type : VDResamplerRegeneratePlots
//     InputRootFile    : "merged_TargetStops1809.root"   # hadd of one source's jobs
//     TreeName         : "VDResamplerGenerator/ttree"
//     trainingPlanFile : "Offline/STMMC/fcl/VDResamplerTrainingPlan.fcl"
//     hitSummaryFile   : "etc.mu2e.STMVDResamplerConfigure_VD116_hitSummary.MDC2025au.001430_00000002.txt"
//     MotherRootFile   : "nts.mu2e.STMVDResamplerConfigure_VD116_hitDump.MDC2025au.001430_00000002.root"
//     MotherTreeName   : "VDResamplerTrainingSetup/ttree"
//   }
//   services.TFileService.fileName : "validation_TargetStops1809.root"
//
// Yongyi Wu, Aug. 2026

#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TFile.h"
#include "TTree.h"

#include "Offline/GeneralUtilities/inc/ParameterSetFromFile.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/STMMC/inc/VDResamplerConfigureCommon.hh"   // resolvePlanEntry (plan lookup)
#include "Offline/STMMC/inc/VDResamplerGenerateCommon.hh"
#include "Offline/STMMC/inc/VDResamplerNameHelper.hh"        // parseSummaryFileName
#include "Offline/STMMC/inc/VDResamplerTrainCommon.hh"       // parse*Basis, validateGeometry
#include "Offline/STMMC/inc/VDResamplerTransforms.hh"
#include "Offline/STMMC/inc/VDResamplerValidationPlots.hh"

namespace mu2e {

  class VDResamplerRegeneratePlots : public art::EDAnalyzer {
    public:
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;

      struct Config {
        fhicl::Atom<std::string> InputRootFile{
          Name("InputRootFile"),
          Comment("Generated-sample ROOT file to replay: one job's output, or the hadd of a "
                  "campaign's jobs FOR A SINGLE SOURCE. Sources are handled one run at a "
                  "time; hadd and run this per source.")};
        fhicl::Atom<std::string> TreeName{
          Name("TreeName"),
          Comment("Path of the generated dump tree inside InputRootFile. This is the "
                  "generator's doROOTDump tree, i.e. '<generatorLabel>/ttree'."),
          "VDResamplerGenerator/ttree"};

        // The plan is the single source of truth for each particle's transform and for the
        // VD geometry, exactly as it is for training and for generation.
        fhicl::Atom<std::string> trainingPlanFile{
          Name("trainingPlanFile"),
          Comment("The SAME training plan fhicl the models were trained with. Supplies each "
                  "particle's SBDMmomentumBasis / SBDMpositionBasis (resolved for this "
                  "source) and the common_training_config geometry, so none of it has to be "
                  "restated here."),
          "Offline/STMMC/fcl/VDResamplerTrainingPlan.fcl"};
        fhicl::Atom<std::string> hitSummaryFile{
          Name("hitSummaryFile"),
          Comment("The hit summary of the source these samples were generated from. Only the "
                  "FILE NAME is parsed (sourceTag/versionTag/runNumber); the contents are not "
                  "read, so a path that no longer resolves is fine as long as the name is "
                  "intact. The source tag it yields is what selects the plan entries.")};

        // The mother side is re-read from the original source file: the dump tree holds only
        // what was generated, so the reference has to come from where it always came from.
        // The summary does not name it, hence this entry.
        fhicl::Atom<std::string> MotherRootFile{
          Name("MotherRootFile"),
          Comment("Source (mother) ROOT file for this source — the SAME file the generation "
                  "job's resamplerSourceRootFiles named for it, otherwise the comparison is "
                  "against a different reference.")};
        fhicl::Atom<std::string> MotherTreeName{
          Name("MotherTreeName"),
          Comment("Tree path within MotherRootFile."),
          "VDResamplerTrainingSetup/ttree"};

        fhicl::OptionalAtom<int> SourceIndex{
          Name("SourceIndex"),
          Comment("Normally unset. The tree's summaryIndex is expected to hold ONE value "
                  "(a single-source file); set this to select one source's entries out of a "
                  "multi-source tree and discard the rest. Only do that when the discarded "
                  "sources are genuinely not wanted — the plots then describe a subset of the "
                  "file, and the POT scaling of the mix no longer applies to them.")};

        // Transform scales. Unlike the geometry these are not part of common_training_config
        // — they are compile-time constants shared by training and generation — so they are
        // exposed here only for the rare case of a non-default transform.
        fhicl::Atom<double> X0{
          Name("X0"), Comment("Detector-axis x origin [mm]"), VDResampler::kX0};
        fhicl::Atom<double> Y0{
          Name("Y0"), Comment("Detector-axis y origin [mm]"), VDResampler::kY0};
        fhicl::Atom<double> T0{
          Name("T0"), Comment("Time transform offset"), VDResampler::kT0};
        fhicl::Atom<double> TScale{
          Name("TScale"), Comment("Time transform scale"), VDResampler::kTScale};
        fhicl::Atom<double> P0{
          Name("P0"), Comment("Momentum transform scale [MeV/c]"), VDResampler::kP0};
      };

      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit VDResamplerRegeneratePlots(const Parameters& conf);
      void analyze(const art::Event&) override {}
      void endJob() override;

    private:
      // The transform for one particle, resolved from the plan for (sourceTag_, pdgId).
      struct SetBasis {
        VDResampler::MomentumBasis momentum = VDResampler::MomentumBasis::V1_CylindricalTransformed;
        VDResampler::PositionBasis position = VDResampler::PositionBasis::V1_Atanh;
      };

      SetBasis basisFor(int pdgId) const;
      std::string labelFor(int pdgId) const;

      std::string inputRootFile_, treeName_, trainingPlanFile_;
      std::string motherFile_, motherTree_;
      fhicl::ParameterSet trainingPlan_;

      // Provenance recovered from the summary file NAME, mirroring what loadSourceSummary
      // recovers in the generator. Fixed for the whole job (one source per run).
      std::string sourceTag_, versionTag_;
      int runNumber_ = 0;

      bool haveSourceIndex_ = false;
      int  sourceIndex_ = 0;

      unsigned long virtualDetectorID_ = 116;
      VDResampler::InverseParams ip_{};
      GlobalConstantsHandle<ParticleDataList> pdt_;
  };

  VDResamplerRegeneratePlots::VDResamplerRegeneratePlots(const Parameters& conf) :
    art::EDAnalyzer(conf),
    inputRootFile_   (conf().InputRootFile()),
    treeName_        (conf().TreeName()),
    trainingPlanFile_(conf().trainingPlanFile()),
    motherFile_      (conf().MotherRootFile()),
    motherTree_      (conf().MotherTreeName())
  {
    if (motherFile_.empty())
      throw cet::exception("VDResamplerRegeneratePlots")
        << "MotherRootFile is empty; there is no reference distribution to compare against.";

    if (trainingPlanFile_.empty())
      throw cet::exception("VDResamplerRegeneratePlots")
        << "trainingPlanFile is required: it is where each particle's momentum and position "
        << "basis come from.";
    trainingPlan_ = ParameterSetFromFile(trainingPlanFile_).pSet();

    // Geometry from the plan's common_training_config, so it cannot drift from what training
    // used. These three define the coordinate system BOTH sides are compared in.
    if (!trainingPlan_.has_key("common_training_config"))
      throw cet::exception("VDResamplerRegeneratePlots")
        << "trainingPlanFile " << trainingPlanFile_ << " has no 'common_training_config' table.";
    const fhicl::ParameterSet common =
      trainingPlan_.get<fhicl::ParameterSet>("common_training_config");
    virtualDetectorID_ = static_cast<unsigned long>(common.get<int>("VirtualDetectorID"));
    ip_ = VDResampler::InverseParams{conf().X0(), conf().Y0(), conf().T0(),
                                     conf().TScale(), conf().P0(),
                                     common.get<double>("VDr"),
                                     common.get<double>("VDz0")};
    VDResampler::validateGeometry(ip_.VDr, ip_.VDz0, "VDResamplerRegeneratePlots");

    // Recover the source's provenance from the summary file NAME, exactly as the generator's
    // loadSourceSummary does — the shared parser is what keeps the two readings identical.
    const std::string summary = conf().hitSummaryFile();
    int summaryVDID = 0;
    if (!VDResampler::parseSummaryFileName(summary, versionTag_, summaryVDID,
                                           sourceTag_, runNumber_))
      throw cet::exception("VDResamplerRegeneratePlots")
        << "Cannot parse hit summary file name \"" << summary << "\". It must follow the "
        << "VDResamplerConfigure convention, which is what carries the source tag, version "
        << "and run number this module needs.";
    if (static_cast<unsigned long>(summaryVDID) != virtualDetectorID_)
      throw cet::exception("VDResamplerRegeneratePlots")
        << "Hit summary \"" << summary << "\" is for VD " << summaryVDID
        << ", but common_training_config.VirtualDetectorID in " << trainingPlanFile_
        << " is " << virtualDetectorID_ << ". The mother selection would then read a "
        << "different detector than the samples were generated on.";

    int index = 0;
    if (conf().SourceIndex(index)) {
      haveSourceIndex_ = true;
      sourceIndex_ = index;
    }
  }

  // This particle's transform, from the plan entry for (sourceTag_, pdgId). resolvePlanEntry
  // applies the plan's documented pdg -> source -> default fallback, so a particle inherits
  // whatever its models were actually trained with.
  VDResamplerRegeneratePlots::SetBasis
  VDResamplerRegeneratePlots::basisFor(int pdgId) const {
    const fhicl::ParameterSet entry = VDResampler::resolvePlanEntry(
        trainingPlan_, sourceTag_, pdgId, "VDResamplerRegeneratePlots");
    SetBasis basis;
    // The defaults match the training modules' own fhicl defaults, so a plan entry that
    // omits a basis key resolves the same way training resolved it.
    basis.momentum = VDResampler::parseMomentumBasis(
        entry.get<std::string>("SBDMmomentumBasis", "V2_PTOT_SLOPES"),
        "VDResamplerRegeneratePlots");
    basis.position = VDResampler::parsePositionBasis(
        entry.get<std::string>("SBDMpositionBasis", "V1_ATANH"),
        "VDResamplerRegeneratePlots");
    return basis;
  }

  // Directory name, built from the same provenance components the generator's label uses
  // (source tag, version, run, pdg) so a regenerated directory is recognisable beside the
  // per-job ones it replaces. Every object of this particle's comparison — histograms,
  // canvases and its validationMetrics tree — is booked inside it, so particles cannot mix.
  std::string VDResamplerRegeneratePlots::labelFor(int pdgId) const {
    return sourceTag_ + "_" + versionTag_
         + "_run" + VDResampler::zeroPad(runNumber_, 6)
         + "_pdg" + VDResampler::pdgFileToken(pdgId);
  }

  void VDResamplerRegeneratePlots::endJob() {
    // TFile::Open (not the TFile ctor) so xroot:// paths work, matching every other reader
    // in this package.
    auto fin = std::unique_ptr<TFile>{TFile::Open(inputRootFile_.c_str(), "READ")};
    if (!fin || fin->IsZombie())
      throw cet::exception("VDResamplerRegeneratePlots")
        << "Cannot open generated ROOT file: " << inputRootFile_;
    TTree* ttree = dynamic_cast<TTree*>(fin->Get(treeName_.c_str()));
    if (!ttree)
      throw cet::exception("VDResamplerRegeneratePlots")
        << "Cannot find TTree '" << treeName_ << "' in " << inputRootFile_
        << ". This must be the generator's doROOTDump tree, e.g. "
        << "'VDResamplerGenerator/ttree'.";

    // Branch schema of the generator's dump tree. z and E are present but not read: z is
    // VDz0 for every sample (no distribution to compare), and E is recoverable from the
    // momentum, so neither is a comparison dimension. generationMode is likewise recorded
    // but does not partition the comparison — a particle's basis is the same either way.
    double x = 0., y = 0., t = 0., px = 0., py = 0., pz = 0.;
    int summaryIndex = 0, pdgId = 0;
    ttree->SetBranchAddress("x",            &x);
    ttree->SetBranchAddress("y",            &y);
    ttree->SetBranchAddress("time",         &t);
    ttree->SetBranchAddress("px",           &px);
    ttree->SetBranchAddress("py",           &py);
    ttree->SetBranchAddress("pz",           &pz);
    ttree->SetBranchAddress("pdgId",        &pdgId);
    ttree->SetBranchAddress("summaryIndex", &summaryIndex);

    const Long64_t nEntries = ttree->GetEntries();
    if (nEntries == 0)
      throw cet::exception("VDResamplerRegeneratePlots")
        << "Tree '" << treeName_ << "' in " << inputRootFile_ << " is empty; there is "
        << "nothing to regenerate.";

    // ---- pass 1: confirm the file really is single-source, and find the particles in it.
    // Booking a set is expensive (it re-reads the whole mother file), so the particles are
    // discovered first and the mother is read once per particle, rather than lazily on first
    // sight while the fill loop is running.
    std::map<int, long> pdgCounts;      // pdgId -> entries kept
    std::set<int> summaryIndices;       // every distinct summaryIndex seen, for the check below
    long skippedOtherSource = 0;
    for (Long64_t i = 0; i < nEntries; ++i) {
      ttree->GetEntry(i);
      summaryIndices.insert(summaryIndex);
      if (haveSourceIndex_ && summaryIndex != sourceIndex_) { ++skippedOtherSource; continue; }
      ++pdgCounts[pdgId];
    }

    // A multi-source tree has no single mother file to be compared against, so it is rejected
    // rather than silently compared against whichever source MotherRootFile happens to name.
    // SourceIndex is the deliberate opt-in for reading one source out of such a file.
    if (!haveSourceIndex_ && summaryIndices.size() > 1) {
      std::ostringstream seen;
      for (int idx : summaryIndices) seen << ' ' << idx;
      throw cet::exception("VDResamplerRegeneratePlots")
        << "Tree '" << treeName_ << "' carries " << summaryIndices.size()
        << " distinct summaryIndex values (" << seen.str() << " ), so it mixes several data "
        << "sources. This module handles ONE source per run: every set in it shares the "
        << "source's mother file, plan entries and geometry. hadd each source's jobs "
        << "separately and run this once per source, or set SourceIndex to select one "
        << "source's entries out of this file and discard the rest.";
    }
    if (haveSourceIndex_ && pdgCounts.empty())
      throw cet::exception("VDResamplerRegeneratePlots")
        << "SourceIndex is " << sourceIndex_ << ", but no entry in '" << treeName_
        << "' has that summaryIndex; nothing would be plotted.";

    // ---- book one plot set per particle, each in its own subdirectory and its own basis.
    // Everything lands under this module's own "validation" directory, so the output file
    // mirrors the generator's layout and existing drawing scripts find the plots in the
    // same place.
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory validationDir = tfs->mkdir("validation");

    std::map<int, std::unique_ptr<VDResampler::ValidationPlots>> plots;
    std::map<int, SetBasis> bases;
    std::ostringstream report;
    report << "Replaying source \"" << sourceTag_ << "\" (" << versionTag_ << ", run "
           << runNumber_ << ") from " << inputRootFile_ << ":" << treeName_ << ": "
           << nEntries << " entr(y/ies), " << pdgCounts.size()
           << " particle(s), each in the basis " << trainingPlanFile_ << " assigns it:";

    for (const auto& [pdg, count] : pdgCounts) {
      const SetBasis basis = basisFor(pdg);
      const std::string label = labelFor(pdg);
      report << "\n  " << label << ": " << count << " sample(s), momentum="
             << VDResampler::momentumBasisName(basis.momentum) << ", position="
             << VDResampler::positionBasisName(basis.position);

      auto set = std::make_unique<VDResampler::ValidationPlots>();
      // No TransformedStatsBySlot is passed: those come from the loaded models' stored
      // training statistics, and this module deliberately loads no models — it replays a
      // tree. The transformed axes therefore take their STATIC fallback ranges, which are
      // wide enough for any sample but coarser over the populated region than the
      // stats-sized axes a live job produces. The physical axes are unaffected (they are
      // always static), so the physical comparison is identical either way.
      //
      // The subdirectory handed in here is what isolates this particle: ValidationPlots
      // books its histograms, its canvases AND its validationMetrics tree into it, so no
      // object of one particle's comparison can land in another's.
      set->book(validationDir.mkdir(label), label, basis.momentum, basis.position,
                pdg, pdt_->particle(pdg).mass(),
                motherFile_, motherTree_, virtualDetectorID_, ip_,
                "VDResamplerRegeneratePlots");
      plots.emplace(pdg, std::move(set));
      bases.emplace(pdg, basis);
    }
    mf::LogInfo("VDResamplerRegeneratePlots") << report.str();

    // ---- pass 2: refill the generated side, one entry at a time.
    VDResampler::PzFallbackStats pzFallback;
    long skippedNonPositivePz = 0;
    for (Long64_t i = 0; i < nEntries; ++i) {
      ttree->GetEntry(i);
      if (haveSourceIndex_ && summaryIndex != sourceIndex_) continue;

      // The forward transform's slope division needs pz > 0, the same condition the mother
      // selection imposes. A generated sample CAN come out with pz <= 0 (nothing bounds the
      // model's output), so those are dropped rather than pushed through the transform
      // against a safety floor — and counted, since a large fraction would itself be a
      // finding about the model.
      if (pz <= 0.0) { ++skippedNonPositivePz; continue; }

      // Reconstruct the transformed coordinates from the physical ones, in THIS particle's
      // basis. This is the round trip described in the header comment: the live job had the
      // model's own transformed output to hand, a replay does not. The mother side goes
      // through this identical call (ValidationPlots::fillMother), so both sides stay in the
      // same coordinates.
      const SetBasis& basis = bases.at(pdgId);
      VDResampler::GeneratedTransformed g;
      g.basis    = basis.momentum;
      g.posBasis = basis.position;
      VDResampler::forwardTransformSample(
          x, y, ip_.VDz0, t, px, py, pz,
          ip_.x0, ip_.y0, ip_.t0, ip_.tScale, ip_.p0, ip_.VDr, ip_.VDz0,
          g.xTrans, g.yTrans, g.tTrans, g.m0, g.m1, g.m2,
          basis.momentum, &pzFallback, basis.position);

      plots.at(pdgId)->fillGenerated(g, x, y, t, px, py, pz);
    }

    fin->Close();

    if (skippedOtherSource > 0)
      mf::LogInfo("VDResamplerRegeneratePlots")
        << "SourceIndex " << sourceIndex_ << " selected: " << skippedOtherSource << " of "
        << nEntries << " entr(y/ies) belong to other sources and were discarded. The plots "
        << "describe only the selected source.";

    if (skippedNonPositivePz > 0)
      mf::LogWarning("VDResamplerRegeneratePlots")
        << skippedNonPositivePz << " generated sample(s) had pz <= 0 and were dropped: the "
        << "forward transform's slope division is undefined there, and the mother selection "
        << "excludes such hits too, so keeping them would compare against a reference that "
        << "cannot contain them. A large fraction here is a statement about the model, not "
        << "about this module.";

    if (pzFallback.count > 0) {
      std::ostringstream oss;
      oss << "pz fell below kPzSafetyEpsilon (" << VDResampler::kPzSafetyEpsilon << ") in "
          << pzFallback.count << " generated sample(s); the floor was used in a slope "
          << "division. First " << pzFallback.firstValues.size() << " offending pz value(s):";
      for (double v : pzFallback.firstValues) oss << ' ' << v;
      mf::LogWarning("VDResamplerRegeneratePlots") << oss.str();
    }

    // ---- finalize: rebin to the MERGED statistics, normalize, score, draw. This is the
    // same call the generator makes at endJob, so the output is the per-job plots' exact
    // counterpart — only with every job's samples behind it.
    for (auto& [pdg, set] : plots) set->finalize();

    mf::LogInfo("VDResamplerRegeneratePlots")
      << "Wrote " << plots.size() << " regenerated comparison set(s) under 'validation/' for "
      << "source \"" << sourceTag_ << "\".";
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::VDResamplerRegeneratePlots)
