// VDResamplerGenerateFromModel_module.cc
// This module takes either a single all-at-once model file or a pair of stage-1/stage-2
// model files and produces new virtual-detector-like GenParticles from the trained
// score-based diffusion model(s).
// Yongyi Wu, Mar. 2026

// stdlib includes
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Offline/MachineLearningTools/inc/ScoreBasedDiffusionModel.hh"

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/STMMC/inc/VDResamplerTransforms.hh"
#include "Offline/STMMC/inc/VDResamplerPtotResampler.hh"
#include "Offline/STMMC/inc/VDResamplerGenerateCommon.hh"
#include "Offline/STMMC/inc/VDResamplerNameHelper.hh"
#include "Offline/STMMC/inc/VDResamplerValidationPlots.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

namespace mu2e {

  namespace {
    // Utility function to extract PDG ID from the model filename, which is expected to contain a substring like "pdg[optional m][digits]", e.g. "pdg13" for muons, "pdgm13" for muon pluses.
    int loadPDGIdFromFileName(const std::string& fileName) {
      const size_t pdgPos = fileName.find("pdg");
      if (pdgPos == std::string::npos) {
        throw cet::exception("VDResamplerGenerateFromModel")
          << "Cannot infer PDG ID from model filename: " << fileName;
      }
      size_t pos = pdgPos + 3;
      bool negative = false;
      if (pos < fileName.size() && (fileName[pos] == 'm' || fileName[pos] == '-')) {
        negative = true;
        ++pos;
      }
      const size_t startDigits = pos;
      while (pos < fileName.size() && std::isdigit(static_cast<unsigned char>(fileName[pos]))) {
        ++pos;
      }
      if (startDigits == pos) {
        throw cet::exception("VDResamplerGenerateFromModel")
          << "Cannot infer PDG ID from model filename: " << fileName;
      }
      const int magnitude = std::stoi(fileName.substr(startDigits, pos - startDigits));
      return negative ? -magnitude : magnitude;
    }

    // Reject a checkpoint whose stored ModelLayout doesn't match the generation slot it is
    // loaded into (an all-at-once file dropped into a stage slot, or vice versa, otherwise
    // generates silently through the wrong inverse transform). The layout tag is only
    // unambiguous for the V2 bases: V1 two-stage stages and the V1 all-at-once model all
    // share layout AllAtOnce6D (the tag predates the stage layouts), so there is nothing to
    // check there and we skip. A pre-tag (v<=6) file reads back as tag 0 == V1 AllAtOnce6D,
    // which also lands here and is correctly left unenforced.
    void requireLayout(const ScoreBasedDiffusionModel& model,
                       VDResampler::ModelLayout expected, const std::string& file,
                       const char* slot) {
      const int tag = model.basisTag();
      if (VDResampler::unpackMomentumBasis(tag) == VDResampler::MomentumBasis::V1_CylindricalTransformed)
        return; // V1 layouts are ambiguous; nothing to enforce
      if (VDResampler::unpackModelLayout(tag) != expected)
        throw cet::exception("VDResamplerGenerateFromModel")
          << "Checkpoint " << file << " loaded into the " << slot << " slot has "
          << VDResampler::basisTagToString(tag) << ", but that slot requires layout "
          << VDResampler::modelLayoutName(expected)
          << ". The model file parameter points at the wrong model.";
    }
  }

  class VDResamplerGenerateFromModel : public art::EDProducer {
    public:
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      struct Config {
        fhicl::Atom<bool> useTwoStageModel{
          Name("useTwoStageModel"),
          Comment("If true, use stage-1/stage-2 models. If false, use the all-at-once 6D model."),
          true
        };
        fhicl::Atom<std::string> stage1ModelFile{
          Name("stage1ModelFile"),
          Comment("Model filename (.dat, or legacy .bin/.csv) for the stage-1 model parameters"),
          ""
        };
        fhicl::Atom<std::string> stage2ModelFile{
          Name("stage2ModelFile"),
          Comment("Model filename (.dat, or legacy .bin/.csv) for the stage-2 model parameters"),
          ""
        };
        fhicl::Atom<std::string> allAtOnceModelFile{
          Name("allAtOnceModelFile"),
          Comment("Model filename (.dat, or legacy .bin/.csv) for the all-at-once 6D model parameters"),
          ""
        };
        fhicl::Atom<std::string> SBDMstage1Method{
          Name("SBDMstage1Method"),
          Comment("V2 two-stage stage-1 pTotal source: DIFFUSION (trained 1-D model) "
                  "or a non-diffusion resampler INVERSE_CDF / SPLINE_CDF / KDE. "
                  "Applies ONLY when useTwoStageModel is true — an all-at-once 6-D model draws "
                  "pTotal within its single sample and has no stage-1 to source, so this is "
                  "ignored (with a warning) in that mode."),
          "DIFFUSION"
        };
        // Peak tagging. Required (and only used) when the stage-2 model was TRAINED with the
        // peak-tag condition dim — its basisTag records that, and generateTwoStage throws if
        // the model expects tags and none are given. These MUST be the same line centers and
        // half-widths, in the same order, that the model was trained with (SBDMpeakTagCenters
        // / SBDMpeakTagHalfWidths in the train module): the network input is the line's INDEX
        // in this list, so a reordered or different list relabels the classes.
        fhicl::Sequence<double> SBDMpeakTagCenters{
          Name("SBDMpeakTagCenters"),
          Comment("Per-tag line center in RAW MeV/c on pTotal; must match the trained model's list."),
          std::vector<double>()
        };
        fhicl::Sequence<double> SBDMpeakTagHalfWidths{
          Name("SBDMpeakTagHalfWidths"),
          Comment("Per-tag inclusive half-window in RAW MeV/c; must match the trained model's list."),
          std::vector<double>()
        };
        fhicl::Atom<std::string> resamplerSourceRootFile{
          Name("resamplerSourceRootFile"),
          Comment("Source ROOT file for the pTotal resampler (required when SBDMstage1Method != DIFFUSION)."),
          ""
        };
        fhicl::Atom<std::string> resamplerSourceTreeName{
          Name("resamplerSourceTreeName"),
          Comment("TTree name in resamplerSourceRootFile."),
          "ttree"
        };
        fhicl::Atom<unsigned long> VirtualDetectorID{
          Name("VirtualDetectorID"),
          Comment("VD id selection for the resampler source (must match training)."),
          116
        };
        fhicl::Atom<bool> useHeun{
          Name("useHeun"),
          Comment("If true, use Heun's method for reverse diffusion. Otherwise use Euler."),
          true
        };
        fhicl::Atom<bool> useSDE{
          Name("useSDE"),
          Comment("If true, use SDE solver. Otherwise use ODE solver."),
          true
        };
        fhicl::Atom<bool> useEMANetworkIfAvailable{
          Name("useEMANetworkIfAvailable"),
          Comment("If true, use the EMA network for inference when available. Pass false to force the base score network."),
          true
        };
        fhicl::Atom<int> diffusionSteps{
          Name("diffusionSteps"),
          Comment("Number of reverse-diffusion steps used for sampling"),
          200
        };
        fhicl::Atom<double> sdeToOdeSigmaThreshold{
          Name("sdeToOdeSigmaThreshold"),
          Comment("Switch from SDE to ODE when sigma falls below this threshold (-1 = always use useSDE setting)"),
          -1.0
        };
        fhicl::Atom<double> VDz0{
          Name("VDz0"),
          Comment("Nominal z coordinate of the virtual detector"),
          37700.39
        };
        fhicl::Atom<double> VDr{
          Name("VDr"),
          Comment("Virtual detector radius used in the training transform"),
          2000.0
        };
        fhicl::Atom<bool> doROOTDump{
          Name("doROOTDump"),
          Comment("Whether to dump generated samples into a ROOT tree"),
          false
        };
        fhicl::Atom<bool> doValidationPlots{
          Name("doValidationPlots"),
          Comment("Whether to write per-dimension generated-vs-mother comparison histograms and "
                  "their W1/JSD/TV/KS metrics into the ROOT file. Requires doROOTDump, and requires "
                  "resamplerSourceRootFile (the mother distribution is read from it) regardless of "
                  "SBDMstage1Method."),
          false
        };
      };

      using Parameters = art::EDProducer::Table<Config>;

      explicit VDResamplerGenerateFromModel(const Parameters& conf);
      void produce(art::Event& event) override;
      void endJob() override;

    private:

      art::RandomNumberGenerator::base_engine_t& engine_;
      CLHEP::RandFlat randFlat_;
      CLHEP::RandGaussQ randGaussQ_;
      bool useTwoStageModel_;
      std::string stage1ModelFile_;
      std::string stage2ModelFile_;
      std::string allAtOnceModelFile_;
      bool useHeun_;
      bool useSDE_;
      bool useEMANetworkIfAvailable_;
      int diffusionSteps_;
      double sdeToOdeSigmaThreshold_;
      double VDz0_;
      double VDr_;
      bool doROOTDump_;
      bool doValidationPlots_;
      GlobalConstantsHandle<ParticleDataList> pdt_;

      std::unique_ptr<ScoreBasedDiffusionModel> allAtOnceModel_;
      std::unique_ptr<ScoreBasedDiffusionModel> stage1Model_;
      std::unique_ptr<ScoreBasedDiffusionModel> stage2Model_;

      // Momentum basis + per-model layout, recovered from the loaded model(s)'
      // opaque basisTag() so the inverse transform auto-selects (no fcl needed).
      VDResampler::MomentumBasis momBasis_ = VDResampler::MomentumBasis::V1_CylindricalTransformed;
      VDResampler::PositionBasis posBasis_ = VDResampler::PositionBasis::V1_Atanh;

      // V2 stage-1 pTotal source. When != DIFFUSION the 1-D stage-1 diffusion model is
      // replaced by ptotResampler_ (built from a required ROOT source file at ctor time).
      VDResampler::Stage1Method stage1Method_ = VDResampler::Stage1Method::DIFFUSION;
      VDResampler::PtotResampler ptotResampler_;

      // Monoenergetic lines whose index becomes the stage-2 class label. Given explicitly here
      // (rather than resolved from a training plan as VDResamplerGenerateMix does) because this
      // module drives one model directly and carries no (pdg, source) plan coordinates. Empty
      // unless the stage-2 model was trained peak-tagged; generateTwoStage checks that against
      // the model's own basisTag.
      std::vector<VDResampler::PeakTag> peakTags_;

      int pdgId_ = 0;

      // Detector-center parameters used in the same transform as training.
      double x0_ = VDResampler::kX0;
      double y0_ = VDResampler::kY0;
      double t0_ = VDResampler::kT0;
      double tScale_ = VDResampler::kTScale;
      double p0_ = VDResampler::kP0;

      // Per-dimension generated-vs-mother comparison, written when doValidationPlots_.
      VDResampler::ValidationPlots validationPlots_;

      // Variables for optional ROOT dump.
      TTree* outTree_ = nullptr;
      double x_gen_ = 0.0;
      double y_gen_ = 0.0;
      double z_gen_ = 0.0;
      double t_gen_ = 0.0;
      double px_gen_ = 0.0;
      double py_gen_ = 0.0;
      double pz_gen_ = 0.0;
      double mass_gen_ = 0.0;
      double E_gen_ = 0.0;
    };

  VDResamplerGenerateFromModel::VDResamplerGenerateFromModel(const Parameters& conf)
    : art::EDProducer{conf},
      engine_(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
      randFlat_(engine_),
      randGaussQ_(engine_),
      useTwoStageModel_(conf().useTwoStageModel()),
      stage1ModelFile_(conf().stage1ModelFile()),
      stage2ModelFile_(conf().stage2ModelFile()),
      allAtOnceModelFile_(conf().allAtOnceModelFile()),
      useHeun_(conf().useHeun()),
      useSDE_(conf().useSDE()),
      useEMANetworkIfAvailable_(conf().useEMANetworkIfAvailable()),
      diffusionSteps_(conf().diffusionSteps()),
      sdeToOdeSigmaThreshold_(conf().sdeToOdeSigmaThreshold()),
      VDz0_(conf().VDz0()),
      VDr_(conf().VDr()),
      doROOTDump_(conf().doROOTDump()),
      doValidationPlots_(conf().doValidationPlots()) {

    produces<GenParticleCollection>();

    if (VDr_ <= 0.0) {
      throw cet::exception("VDResamplerGenerateFromModel")
        << "VDr must be positive (got " << VDr_ << "); "
        << "rho = r/VDr would produce inf/NaN in generated samples.";
    }
    if (!std::isfinite(VDz0_)) {
      throw cet::exception("VDResamplerGenerateFromModel")
        << "VDz0 must be finite (got " << VDz0_ << ").";
    }

    stage1Method_ = VDResampler::parseStage1Method(conf().SBDMstage1Method(), "VDResamplerGenerateFromModel");
    const bool useResampler = (stage1Method_ != VDResampler::Stage1Method::DIFFUSION);

    // Peak-tag lines. Validated to the same rules as the train module's assemblePeakTags
    // (equal lengths, positive values, disjoint windows) so a typo is caught at setup rather
    // than becoming a silently wrong class label on every generated event. Whether they are
    // actually REQUIRED is decided by the stage-2 model's own basisTag, in generateTwoStage.
    {
      const std::vector<double> centers    = conf().SBDMpeakTagCenters();
      const std::vector<double> halfWidths = conf().SBDMpeakTagHalfWidths();
      if (centers.size() != halfWidths.size())
        throw cet::exception("VDResamplerGenerateFromModel")
          << "SBDMpeakTagCenters (" << centers.size() << " entries) and SBDMpeakTagHalfWidths ("
          << halfWidths.size() << ") must have the same length.";
      for (size_t k = 0; k < centers.size(); ++k) {
        if (!(centers[k] > 0.0) || !(halfWidths[k] > 0.0))
          throw cet::exception("VDResamplerGenerateFromModel")
            << "Peak tag " << k << " has center " << centers[k] << " and half-width "
            << halfWidths[k] << "; both must be > 0 (MeV/c).";
        peakTags_.push_back(VDResampler::PeakTag{centers[k], halfWidths[k]});
      }
      for (size_t a = 0; a < peakTags_.size(); ++a)
        for (size_t b = a + 1; b < peakTags_.size(); ++b)
          if (std::abs(peakTags_[a].center - peakTags_[b].center)
              <= peakTags_[a].halfWidth + peakTags_[b].halfWidth)
            throw cet::exception("VDResamplerGenerateFromModel")
              << "Peak tag windows " << a << " and " << b << " overlap; they must be disjoint "
              << "so the class label does not depend on the order they are configured in.";
    }

    if (useTwoStageModel_) {
      // stage1 may be supplied either by a trained 1-D diffusion model OR by the
      // pTotal resampler. Stage2 is always required.
      if (stage2ModelFile_.empty()) {
        throw cet::exception("VDResamplerGenerateFromModel")
          << "Two-stage generation requires stage2ModelFile.";
      }
      if (!useResampler && stage1ModelFile_.empty()) {
        throw cet::exception("VDResamplerGenerateFromModel")
          << "Two-stage generation with SBDMstage1Method=DIFFUSION requires stage1ModelFile.";
      }

      stage2Model_ = std::make_unique<ScoreBasedDiffusionModel>(
        ScoreBasedDiffusionModel::loadModel(randFlat_, randGaussQ_, stage2ModelFile_)
      );
      momBasis_ = VDResampler::unpackMomentumBasis(stage2Model_->basisTag());
      posBasis_ = VDResampler::unpackPositionBasis(stage2Model_->basisTag());
      mf::LogInfo("VDResamplerGenerateFromModel")
        << "Loaded stage-2 model " << stage2ModelFile_ << ": "
        << VDResampler::basisTagToString(stage2Model_->basisTag());
      requireLayout(*stage2Model_, VDResampler::ModelLayout::TwoStageStage2_5D,
                    stage2ModelFile_, "stage-2");
      pdgId_ = loadPDGIdFromFileName(stage2ModelFile_);

      // The peak label is only as sharp as the pTotal it is derived from. A DIFFUSION stage-1
      // is a continuous model over log(pTotal/p0): it cannot reproduce a line whose width is a
      // keV, so it spreads those events out and the tag fires on a smeared population rather
      // than on the line — the tagged fraction will not match the truth.
      //
      // A QUALITY warning, not an error: the stage-2 model is trained and used correctly either
      // way (its training labels came from true hits), and a diffusion stage-1 is a legitimate
      // choice. But peak tagging exists precisely to preserve a narrow line, so pairing it with
      // the one stage-1 method that cannot is almost always unintended. (The complementary
      // error — a tagged model with no tags configured — is raised by generateTwoStage.)
      if (VDResampler::unpackPeakTagged(stage2Model_->basisTag()) && !useResampler)
        mf::LogWarning("VDResamplerGenerateFromModel")
          << "Stage-2 model " << stage2ModelFile_ << " was trained with peak tagging, but "
          << "SBDMstage1Method=DIFFUSION. The stage-1 diffusion model smooths over the narrow "
          << "line(s) the tags name, so the fraction of events tagged will not match the source. "
          << "Use a pTotal resampler (INVERSE_CDF / SPLINE_CDF), which draws from the empirical "
          << "distribution and preserves the line exactly.";

      if (useResampler) {
        // Build the pTotal resampler from the required source ROOT file. The basis is
        // read from stage2 (stage1 has no model). The resampler source selection uses
        // the fcl VD id and the pdgId derived from the stage2 file name (so source and
        // model necessarily agree on the particle).
        const std::string srcFile = conf().resamplerSourceRootFile();
        if (srcFile.empty())
          throw cet::exception("VDResamplerGenerateFromModel")
            << "SBDMstage1Method=" << conf().SBDMstage1Method()
            << " requires resamplerSourceRootFile (the pTotal source).";
        ptotResampler_.buildFromRoot(srcFile, conf().resamplerSourceTreeName(),
                                     conf().VirtualDetectorID(), pdgId_,
                                     stage1Method_, "VDResamplerGenerateFromModel");
      } else {
        stage1Model_ = std::make_unique<ScoreBasedDiffusionModel>(
          ScoreBasedDiffusionModel::loadModel(randFlat_, randGaussQ_, stage1ModelFile_)
        );
        mf::LogInfo("VDResamplerGenerateFromModel")
          << "Loaded stage-1 model " << stage1ModelFile_ << ": "
          << VDResampler::basisTagToString(stage1Model_->basisTag());
        requireLayout(*stage1Model_, VDResampler::ModelLayout::TwoStageStage1Ptot1D,
                      stage1ModelFile_, "stage-1");
        // Two loaded models must carry the same basis tag, else the inverse is ambiguous.
        const auto stage1Basis = VDResampler::unpackMomentumBasis(stage1Model_->basisTag());
        if (stage1Basis != momBasis_)
          throw cet::exception("VDResamplerGenerateFromModel")
            << "Stage-1 and stage-2 models disagree on momentum basis (tags "
            << stage1Model_->basisTag() << " vs " << stage2Model_->basisTag() << ").";
      }
    } else {
      if (allAtOnceModelFile_.empty()) {
        throw cet::exception("VDResamplerGenerateFromModel")
          << "All-at-once generation requires allAtOnceModelFile.";
      }
      // stage1Method_ was parsed above (so a typo is still caught) but nothing in this branch
      // consults it: a 6-D all-at-once model draws pTotal as part of the single sample, so it
      // has no stage-1 to source. Setting it to a resampler method here would otherwise be
      // accepted in silence and produce a plain all-at-once job — say so instead.
      if (useResampler) {
        mf::LogWarning("VDResamplerGenerateFromModel")
          << "SBDMstage1Method=" << conf().SBDMstage1Method() << " is IGNORED when "
          << "useTwoStageModel is false: an all-at-once 6-D model has no separate stage-1 "
          << "pTotal step to resample. No pTotal resampler is built, and resamplerSourceRootFile "
          << "is used only if doValidationPlots is on (for the mother distribution). Set "
          << "useTwoStageModel: true if you meant to use the configuration.";
      }

      allAtOnceModel_ = std::make_unique<ScoreBasedDiffusionModel>(
        ScoreBasedDiffusionModel::loadModel(randFlat_, randGaussQ_, allAtOnceModelFile_)
      );
      momBasis_ = VDResampler::unpackMomentumBasis(allAtOnceModel_->basisTag());
      posBasis_ = VDResampler::unpackPositionBasis(allAtOnceModel_->basisTag());
      mf::LogInfo("VDResamplerGenerateFromModel")
        << "Loaded all-at-once model " << allAtOnceModelFile_ << ": "
        << VDResampler::basisTagToString(allAtOnceModel_->basisTag());
      requireLayout(*allAtOnceModel_, VDResampler::ModelLayout::AllAtOnce6D,
                    allAtOnceModelFile_, "all-at-once");
      pdgId_ = loadPDGIdFromFileName(allAtOnceModelFile_);
    }

    z_gen_ = VDz0_;

    if (doROOTDump_) {
      art::ServiceHandle<art::TFileService> tfs;
      outTree_ = tfs->make<TTree>("ttree", "Generated samples from VD resampler model");
      outTree_->Branch("pdgId", &pdgId_, "pdgId/I");
      outTree_->Branch("x", &x_gen_, "x/D");
      outTree_->Branch("y", &y_gen_, "y/D");
      outTree_->Branch("z", &z_gen_, "z/D");
      outTree_->Branch("time", &t_gen_, "time/D");
      outTree_->Branch("px", &px_gen_, "px/D");
      outTree_->Branch("py", &py_gen_, "py/D");
      outTree_->Branch("pz", &pz_gen_, "pz/D");
      outTree_->Branch("E", &E_gen_, "E/D");

      // Validation plots live in this same ROOT file, so they are booked here, under the
      // dump, reusing its TFileService handle.
      if (doValidationPlots_) {
        // The mother distribution comes from the source ROOT file. For a DIFFUSION stage-1
        // job nothing else needs that file, so it is optional there — but validation cannot
        // run without a reference, hence it is mandatory whenever the plots are requested.
        const std::string srcFile = conf().resamplerSourceRootFile();
        if (srcFile.empty())
          throw cet::exception("VDResamplerGenerateFromModel")
            << "doValidationPlots requires resamplerSourceRootFile: it supplies the mother "
            << "distribution the generated samples are compared against.";

        // momBasis_ and pdgId_ are settled above from the loaded model(s), so the plot set is
        // booked for exactly the basis that will be generated.
        const VDResampler::InverseParams ip{x0_, y0_, t0_, tScale_, p0_, VDr_, VDz0_};
        // The models' own training statistics size the transformed axes to where the
        // population actually is, rather than to the static catch-all ranges.
        const VDResampler::TransformedStatsBySlot stats =
          VDResampler::collectTransformedStats(allAtOnceModel_.get(), stage1Model_.get(),
                                               stage2Model_.get(), momBasis_);
        validationPlots_.book(
          tfs->mkdir("validation"),
          "pdg" + VDResampler::pdgFileToken(pdgId_), momBasis_, posBasis_, pdgId_,
          pdt_->particle(pdgId_).mass(),
          srcFile, conf().resamplerSourceTreeName(), conf().VirtualDetectorID(), ip,
          "VDResamplerGenerateFromModel", &stats);
      }
    } else if (doValidationPlots_) {
      throw cet::exception("VDResamplerGenerateFromModel")
        << "doValidationPlots requires doROOTDump (the comparison histograms are written "
        << "into the ROOT dump file).";
    }
  }

  void VDResamplerGenerateFromModel::endJob() {
    validationPlots_.finalize();
  }

  void VDResamplerGenerateFromModel::produce(art::Event& event) {
    auto output = std::make_unique<GenParticleCollection>();

    const VDResampler::SamplerSettings settings{
      useEMANetworkIfAvailable_, useHeun_, useSDE_, diffusionSteps_, sdeToOdeSigmaThreshold_};

    // Generate one transformed sample (shared with VDResamplerGenerateMix).
    VDResampler::GeneratedTransformed g;
    if (useTwoStageModel_) {
      g = VDResampler::generateTwoStage(
        stage1Model_.get(), *stage2Model_, stage1Method_, ptotResampler_,
        randFlat_, randGaussQ_, settings, p0_, "VDResamplerGenerateFromModel", peakTags_);
    } else {
      g = VDResampler::generateAllAtOnce(*allAtOnceModel_, settings, "VDResamplerGenerateFromModel");
    }

    const VDResampler::InverseParams ip{x0_, y0_, t0_, tScale_, p0_, VDr_, VDz0_};
    VDResampler::invertGenerated(g, ip, x_gen_, y_gen_, z_gen_, t_gen_, px_gen_, py_gen_, pz_gen_);

    mass_gen_ = pdt_->particle(pdgId_).mass();
    const CLHEP::Hep3Vector momParticle(px_gen_, py_gen_, pz_gen_);
    const CLHEP::Hep3Vector posParticle(x_gen_, y_gen_, z_gen_);
    const double eTotal = std::sqrt(momParticle.mag2() + mass_gen_ * mass_gen_);
    E_gen_ = eTotal - mass_gen_;
    const CLHEP::HepLorentzVector fourMomParticle(momParticle, eTotal);

    output->emplace_back(
      PDGCode::type(pdgId_),
      GenId::STMDownStreamGenTool,
      posParticle,
      fourMomParticle,
      t_gen_
    );

    event.put(std::move(output));

    if (doROOTDump_) {
      outTree_->Fill();
    }
    // booked() is true only when validation was requested, so it is the standing record of
    // that decision — no separate flag check is needed.
    if (validationPlots_.booked()) {
      // g holds the model's own transformed output; the *_gen_ values are its inversion.
      validationPlots_.fillGenerated(g, x_gen_, y_gen_, t_gen_, px_gen_, py_gen_, pz_gen_);
    }
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::VDResamplerGenerateFromModel)
