// VDResamplerTrainFromRoot_module.cc
// A copy of VDResamplerTrain_module.cc but takes ROOT file input instead of
// StepPointMCCollection and SimParticleCollection from art::Event.
// Using ROOT file input from VDResamplerConfigure_module.cc,
// and saving models at specified epochs as well as the final epoch.
// Yongyi Wu, May 2026

#include <memory>
#include <string>
#include <vector>

#include "Offline/STMMC/inc/VDResamplerTrainCommon.hh"
#include "Offline/STMMC/inc/VDResamplerTransforms.hh"
#include "Offline/SeedService/inc/SeedService.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TFile.h"
#include "TTree.h"

namespace mu2e {
  class VDResamplerTrainFromRoot : public art::EDAnalyzer {
    public:
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      struct Config {
        fhicl::Atom<std::string> InputRootFile{ Name("InputRootFile"), Comment("Input ROOT file with TTree") };
        fhicl::Atom<std::string> TreeName{      Name("TreeName"),      Comment("Name of the TTree in the ROOT file"), "VDResamplerTrainingSetup/ttree" };
        fhicl::Atom<bool>        SBDMuseTwoStageTraining{    Name("SBDMuseTwoStageTraining"),    Comment("If true, train two-stage model; if false, train all 6 dimensions at once."), true };
        fhicl::Atom<std::string> SBDMallAtOnceModelFile{    Name("SBDMallAtOnceModelFile"),    Comment("Model output filename (.dat) for the all-at-once 6D model"),  "" };
        fhicl::Atom<std::string> SBDMstage1ModelFile{        Name("SBDMstage1ModelFile"),        Comment("Model output filename (.dat) for the stage-1 model"),         "" };
        fhicl::Atom<std::string> SBDMstage2ModelFile{        Name("SBDMstage2ModelFile"),        Comment("Model output filename (.dat) for the stage-2 model"),         "" };
        fhicl::Atom<std::string> SBDMloadCheckPointAllAtOnceModelFile{ Name("SBDMloadCheckPointAllAtOnceModelFile"), Comment("Checkpoint file to load for the all-at-once model (.dat, or legacy .bin/.csv)"), "" };
        fhicl::Atom<std::string> SBDMloadCheckPointStage1ModelFile{    Name("SBDMloadCheckPointStage1ModelFile"),    Comment("Checkpoint file to load for the stage-1 model (.dat, or legacy .bin/.csv)"),    "" };
        fhicl::Atom<std::string> SBDMloadCheckPointStage2ModelFile{    Name("SBDMloadCheckPointStage2ModelFile"),    Comment("Checkpoint file to load for the stage-2 model (.dat, or legacy .bin/.csv)"),    "" };
        fhicl::Atom<bool>       SBDMpromoteEMA{                           Name("SBDMpromoteEMA"),                           Comment("Promote EMA weights to network (and reset optimizer) once at the start of training"), false };
        fhicl::Atom<int>    VirtualDetectorID{ Name("VirtualDetectorID"), Comment("Virtual detector ID to select"),    116 };
        fhicl::Atom<double> VDz0{              Name("VDz0"),              Comment("z coordinate of the virtual detector"),      37700.39 };
        fhicl::Atom<double> VDr{               Name("VDr"),               Comment("VD radius"),                                 2000.0 };
        fhicl::Atom<int>    pdgID{             Name("pdgID"),             Comment("PDG ID to select"),                          22 };
        fhicl::Atom<std::string> SBDMmomentumBasis{ Name("SBDMmomentumBasis"), Comment("Momentum transform basis: V1_CYLINDRICAL, V2_PTOT_SLOPES, V2_PTOT_SLOPES_ASINH, V3_PTOT_SLOPES_ASINH_TIME_ASINH"), "V2_PTOT_SLOPES" };
        fhicl::Atom<int>    SBDMtimeEmbeddingDim{ Name("SBDMtimeEmbeddingDim"), Comment("Time embedding dimension"),            0 };
        fhicl::Sequence<int> SBDMinputEmbeddingDims{     Name("SBDMinputEmbeddingDims"),     Comment("Per-state-dim Fourier depth: [] none, [k] broadcast, or length-dim list; each 0 or even >= 2"),     std::vector<int>() };
        fhicl::Sequence<int> SBDMconditionEmbeddingDims{ Name("SBDMconditionEmbeddingDims"), Comment("Per-condition-dim Fourier depth: [] none, [k] broadcast, or length-condDim list; each 0 or even >= 2"), std::vector<int>() };
        fhicl::Atom<int>    SBDMhidden{        Name("SBDMhidden"),        Comment("Hidden layer size"),                         128 };
        fhicl::Atom<int>    SBDMlayers{        Name("SBDMlayers"),        Comment("Number of layers"),                          4 };
        fhicl::Atom<std::string> SBDMoptimizer{ Name("SBDMoptimizer"),   Comment("Optimizer (SGD/ADAM)"),                      "ADAM" };
        fhicl::Atom<double> SBDMadamBeta1{     Name("SBDMadamBeta1"),     Comment("Adam beta1"),                                0.9 };
        fhicl::Atom<double> SBDMadamBeta2{     Name("SBDMadamBeta2"),     Comment("Adam beta2"),                                0.999 };
        fhicl::Atom<double> SBDMadamEps{       Name("SBDMadamEps"),       Comment("Adam epsilon"),                              1e-8 };
        fhicl::Atom<std::string> SBDMnoiseSchedule{ Name("SBDMnoiseSchedule"), Comment("Noise schedule (LINEAR/COSINE/LOGSIG)"), "COSINE" };
        fhicl::Atom<double> SBDMbetaMin{       Name("SBDMbetaMin"),       Comment("Min beta (LINEAR schedule)"),                1e-4 };
        // Continuous VP-SDE: beta(t) is a rate integrated over t in [0,1], not a DDPM per-step beta.
        // The integral ~0.5*betaMax must be O(a few) to fully noise the data by t=1; the old 0.02
        // (DDPM per-step value, meant to sum over ~1000 steps) left sigma(1)~0.1, not ~1.
        fhicl::Atom<double> SBDMbetaMax{       Name("SBDMbetaMax"),       Comment("Max beta (LINEAR schedule)"),                20.0 };
        fhicl::Atom<double> SBDMcosineOffset{  Name("SBDMcosineOffset"),  Comment("Cosine offset"),                             0.008 };
        fhicl::Atom<double> SBDMlogSigMin{     Name("SBDMlogSigMin"),     Comment("Min sigma (LOGSIG schedule)"),               1e-5 };
        fhicl::Atom<double> SBDMlogSigMax{     Name("SBDMlogSigMax"),     Comment("Max sigma (LOGSIG schedule)"),               1.0 };
        fhicl::Atom<std::string> SBDMpredictionTarget{ Name("SBDMpredictionTarget"), Comment("Network regression target: SCORE, EPS, or V (v-prediction)"), "SCORE" };
        // Deprecated: superseded by SBDMpredictionTarget. If set, true->EPS, false->SCORE.
        fhicl::OptionalAtom<bool> SBDMepsPrediction{ Name("SBDMepsPrediction"), Comment("DEPRECATED: use SBDMpredictionTarget. true->EPS, false->SCORE") };
        fhicl::Atom<double> SBDMlossWeightPower{ Name("SBDMlossWeightPower"), Comment("Loss weight power"),                     2.0 };
        fhicl::Atom<int>    SBDMbatchSize{     Name("SBDMbatchSize"),     Comment("Batch size"),                                32 };
        fhicl::Atom<double> SBDMgradientClip{  Name("SBDMgradientClip"),  Comment("Gradient clip threshold"),                   1.0 };
        fhicl::Atom<double> SBDMlearningRate{  Name("SBDMlearningRate"),  Comment("Learning rate"),                             1e-3 };
        fhicl::Atom<bool>   SBDMbiasLowSigma{  Name("SBDMbiasLowSigma"), Comment("Bias training samples towards low-sigma"),   false };
        fhicl::Atom<double> SBDMtLowBound{     Name("SBDMtLowBound"),     Comment("Lower bound clamp on t for training"),       0.0 };
        fhicl::Atom<double> SBDMtFocusLow{      Name("SBDMtFocusLow"),      Comment("Low edge of t focus window"),                                  0.0 };
        fhicl::Atom<double> SBDMtFocusHigh{     Name("SBDMtFocusHigh"),     Comment("High edge of t focus window"),                                 0.0 };
        fhicl::Atom<double> SBDMtFocusFraction{ Name("SBDMtFocusFraction"), Comment("Realized fraction of training samples inside the t focus window; the rest are drawn from its complement (0 = disabled)"), 0.0 };
        fhicl::Atom<bool>   SBDMuseDimWeightController{      Name("SBDMuseDimWeightController"),      Comment("Use adaptive per-dimension gradient weighting"), false };
        fhicl::Atom<double> SBDMdimWeightControllerEMADecay{ Name("SBDMdimWeightControllerEMADecay"), Comment("EMA decay for dimension weight controller"),     0.99 };
        fhicl::Atom<bool>   SBDMuseEMANetwork{      Name("SBDMuseEMANetwork"),      Comment("Maintain EMA copy of network for inference"), false };
        fhicl::Atom<double> SBDMemaNetworkDecay{    Name("SBDMemaNetworkDecay"),    Comment("Decay for EMA network"),                      0.9999 };
        fhicl::Atom<int>    SBDMdiffusionSteps{     Name("SBDMdiffusionSteps"),     Comment("Diffusion steps"),                            200 };
        fhicl::Atom<int>    SBDMtrainingSize{       Name("SBDMtrainingSize"),       Comment("Training data size (-1 = all)"),               -1 };
        // Samples DRAWN per epoch (0 = one full pass; >N cycles the data with reshuffle-on-wrap so an
        // epoch is a fixed amount of optimization work). SBDMtrainingSubsetSizePerEpoch is the
        // deprecated former name, still accepted; the new key wins when set to a non-zero value.
        fhicl::Atom<int>    SBDMsamplesDrawnPerEpoch{ Name("SBDMsamplesDrawnPerEpoch"), Comment("Samples drawn per epoch (0 = full set; >N cycles the data)"), 0 };
        fhicl::Atom<int>    SBDMtrainingSubsetSizePerEpoch{ Name("SBDMtrainingSubsetSizePerEpoch"), Comment("[deprecated: use SBDMsamplesDrawnPerEpoch] Subset size per epoch (0 = full set)"), 0 };
        fhicl::Atom<int>    SBDMtrainingEpochs{     Name("SBDMtrainingEpochs"),     Comment("Number of training epochs"),                  10 };
        fhicl::Sequence<int> SaveEpochs{ Name("SaveEpochs"), Comment("Epochs at which to save checkpoint models"), std::vector<int>() };
        fhicl::Atom<bool>   SBDMautoCurriculumPlanner{ Name("SBDMautoCurriculumPlanner"), Comment("Train each curriculum phase until loss converges instead of for a fixed epoch count"), false };
        fhicl::Atom<int>    SBDMplannerSmoothWindow{    Name("SBDMplannerSmoothWindow"),    Comment("Planner: trailing moving-average window for smoothing per-epoch loss"), 10 };
        fhicl::Atom<double> SBDMplannerMinDelta{        Name("SBDMplannerMinDelta"),        Comment("Planner: relative improvement of smoothed loss counted as progress"), 0.005 };
        fhicl::Atom<int>    SBDMplannerPatience{        Name("SBDMplannerPatience"),        Comment("Planner: consecutive non-improving epochs to declare convergence"), 20 };
        fhicl::Atom<int>    SBDMplannerMinEpochsPerPhase{ Name("SBDMplannerMinEpochsPerPhase"), Comment("Planner: minimum epochs per phase before convergence can fire (auto-raised to >= smoothWindow+patience)"), 30 };
        fhicl::Atom<int>    SBDMplannerMaxEpochsPerPhase{ Name("SBDMplannerMaxEpochsPerPhase"), Comment("Planner: hard cap per phase; hitting it aborts training with a warning"), 100 };
        fhicl::Sequence<int>    SBDMtrainingCurriculumEpochs{             Name("SBDMtrainingCurriculumEpochs"),             Comment("Epochs per curriculum phase"),                std::vector<int>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumLossWeightPower{    Name("SBDMtrainingCurriculumLossWeightPower"),    Comment("Loss weight power per curriculum phase"),     std::vector<double>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumGradientClip{       Name("SBDMtrainingCurriculumGradientClip"),       Comment("Gradient clip per curriculum phase"),         std::vector<double>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumLearningRate{       Name("SBDMtrainingCurriculumLearningRate"),       Comment("Learning rate per curriculum phase"),         std::vector<double>() };
        fhicl::Sequence<bool>   SBDMtrainingCurriculumBiasLowSigma{       Name("SBDMtrainingCurriculumBiasLowSigma"),       Comment("BiasLowSigma per curriculum phase"),          std::vector<bool>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumTLowBound{          Name("SBDMtrainingCurriculumTLowBound"),          Comment("tLowBound per curriculum phase"),             std::vector<double>() };
        fhicl::Sequence<int>    SBDMtrainingCurriculumBatchSize{          Name("SBDMtrainingCurriculumBatchSize"),          Comment("Batch size per curriculum phase"),            std::vector<int>() };
        fhicl::Sequence<bool>   SBDMtrainingCurriculumPromoteEMA{         Name("SBDMtrainingCurriculumPromoteEMA"),         Comment("If true, promote EMA to network when entering this curriculum phase"),               std::vector<bool>() };
        fhicl::Sequence<bool>   SBDMtrainingCurriculumUseDimWeightController{ Name("SBDMtrainingCurriculumUseDimWeightController"), Comment("Enable per-dimension gradient weight controller per curriculum phase"), std::vector<bool>() };
        fhicl::Sequence<bool>   SBDMtrainingCurriculumUsePeakWindowLoss{ Name("SBDMtrainingCurriculumUsePeakWindowLoss"), Comment("Auto-planner: plateau on the peak-window (feature-region) loss instead of the aggregate loss, per curriculum phase (requires SBDMpeakWindow* configured)"), std::vector<bool>() };
        fhicl::Sequence<bool>   SBDMtrainingCurriculumUsePeakSampling{ Name("SBDMtrainingCurriculumUsePeakSampling"), Comment("Per curriculum phase: enable peak-window importance (emphasis) sampling. Empty defaults to true for all phases. When false, the phase draws uniformly (no oversampling, wIS=1) even if SBDMpeakWindow* is configured"), std::vector<bool>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumPeakAlpha{ Name("SBDMtrainingCurriculumPeakAlpha"), Comment("Per-phase override of peak importance-sampling alpha for ALL windows (1=unbiased, <1=up-weight feature); negative = use base SBDMpeakAlphas. Empty = no override."), std::vector<double>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumTFocusLow{          Name("SBDMtrainingCurriculumTFocusLow"),          Comment("t focus window low edge per curriculum phase"),       std::vector<double>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumTFocusHigh{         Name("SBDMtrainingCurriculumTFocusHigh"),         Comment("t focus window high edge per curriculum phase"),      std::vector<double>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumTFocusFraction{     Name("SBDMtrainingCurriculumTFocusFraction"),     Comment("Realized in-window t focus fraction per curriculum phase (rest drawn from complement)"),       std::vector<double>() };
        // Per-phase samples drawn per epoch. SBDMtrainingCurriculumSamplesDrawnPerEpoch is the
        // preferred name; SBDMtrainingCurriculumSubsetSizePerEpoch is the deprecated former name,
        // still accepted. The new key wins when it is non-empty.
        fhicl::Sequence<int>    SBDMtrainingCurriculumSamplesDrawnPerEpoch{ Name("SBDMtrainingCurriculumSamplesDrawnPerEpoch"), Comment("Samples drawn per epoch per curriculum phase (0 = full set; >N cycles the data)"), std::vector<int>() };
        fhicl::Sequence<int>    SBDMtrainingCurriculumSubsetSizePerEpoch{ Name("SBDMtrainingCurriculumSubsetSizePerEpoch"), Comment("[deprecated: use SBDMtrainingCurriculumSamplesDrawnPerEpoch] Training subset size per epoch per curriculum phase (0 = full set)"), std::vector<int>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumMinDelta{ Name("SBDMtrainingCurriculumMinDelta"), Comment("Auto-planner minDelta per curriculum phase (relative smoothed-loss improvement threshold; empty = use SBDMplannerMinDelta)"), std::vector<double>() };
        fhicl::Sequence<double> SBDMdenoiseDiagnosticTs{ Name("SBDMdenoiseDiagnosticTs"), Comment("If non-empty, run the one-step denoising diagnostic at these t values INSTEAD of training"), std::vector<double>() };
        fhicl::Atom<int>        SBDMdenoiseDiagnosticSamples{ Name("SBDMdenoiseDiagnosticSamples"), Comment("Samples per t value for the denoising diagnostic"), 100000 };
        fhicl::Atom<bool>       SBDMdenoiseDiagnosticUseEMA{ Name("SBDMdenoiseDiagnosticUseEMA"), Comment("Network used by the diagnostics: true = EMA network when available, false = base network. Match the generation config's useEMANetworkIfAvailable"), true };
        fhicl::Sequence<double> SBDMpartialReverseT0s{ Name("SBDMpartialReverseT0s"), Comment("If non-empty, run the partial-reverse sampling diagnostic from these t0 values INSTEAD of training"), std::vector<double>() };
        fhicl::Atom<int>        SBDMpartialReverseSamples{ Name("SBDMpartialReverseSamples"), Comment("Samples per t0 value for the partial-reverse diagnostic"), 20000 };
        fhicl::Atom<bool>       SBDMpartialReverseUseHeun{ Name("SBDMpartialReverseUseHeun"), Comment("Use Heun's method in the partial-reverse diagnostic sampler"), true };
        fhicl::Atom<bool>       SBDMpartialReverseUseSDE{ Name("SBDMpartialReverseUseSDE"), Comment("Use the SDE solver in the partial-reverse diagnostic sampler"), true };
        fhicl::Atom<int>        SBDMpartialReverseDiffusionSteps{ Name("SBDMpartialReverseDiffusionSteps"), Comment("Diffusion steps for the partial-reverse diagnostic sampler (-1 = model's configured value)"), -1 };
        fhicl::Atom<double>     SBDMpartialReverseSdeToOdeSigmaThreshold{ Name("SBDMpartialReverseSdeToOdeSigmaThreshold"), Comment("Switch from SDE to ODE when sigma falls below this threshold in the partial-reverse diagnostic sampler (-1 = always use SBDMpartialReverseUseSDE setting)"), -1.0 };
        fhicl::Atom<bool>       SBDMcondLossDiagnostic{ Name("SBDMcondLossDiagnostic"), Comment("If true, write the conditional eps-loss-vs-sigma diagnostic (split by SBDMpeakWindow*) INSTEAD of training"), false };
        fhicl::Atom<int>        SBDMcondLossDiagnosticSamples{ Name("SBDMcondLossDiagnosticSamples"), Comment("Samples for the conditional-loss diagnostic"), 200000 };
        fhicl::Atom<bool>       SBDMfeatureBlockDiagnostic{ Name("SBDMfeatureBlockDiagnostic"), Comment("If true, write the first-layer feature-block weight/grad magnitude diagnostic INSTEAD of training"), false };
        fhicl::Atom<int>        SBDMfeatureBlockDiagnosticSamples{ Name("SBDMfeatureBlockDiagnosticSamples"), Comment("Draws used to accumulate the gradient for the feature-block diagnostic"), 50000 };
        // Peak importance sampling: parallel sequences, one entry per window (all same length; empty = disabled).
        fhicl::Sequence<int>    SBDMpeakWindowDims{  Name("SBDMpeakWindowDims"),  Comment("Per-window normalized state-coordinate index (e.g. pz = 5 all-at-once, 2 stage-2)"), std::vector<int>() };
        fhicl::Sequence<double> SBDMpeakWindowLows{  Name("SBDMpeakWindowLows"),  Comment("Per-window lower edge in TRANSFORMED (pre-z-score) units, e.g. log(pz/p0), inclusive"), std::vector<double>() };
        fhicl::Sequence<double> SBDMpeakWindowHighs{ Name("SBDMpeakWindowHighs"), Comment("Per-window upper edge in TRANSFORMED (pre-z-score) units, exclusive"), std::vector<double>() };
        fhicl::Sequence<double> SBDMpeakGMaxes{      Name("SBDMpeakGMaxes"),      Comment("Per-window sampling-fraction ceiling at low sigma (0<gMax<1; sum<1)"), std::vector<double>() };
        fhicl::Sequence<double> SBDMpeakSigma0s{     Name("SBDMpeakSigma0s"),     Comment("Per-window Gaussian sigma-taper scale (~ feature width)"), std::vector<double>() };
        fhicl::Sequence<double> SBDMpeakAlphas{      Name("SBDMpeakAlphas"),      Comment("Per-window 1=unbiased, <1=up-weight"), std::vector<double>() };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit VDResamplerTrainFromRoot(const Parameters& conf);
      void analyze(const art::Event&) override {}
      void endJob() override;

    private:
      art::RandomNumberGenerator::base_engine_t& engine_;
      CLHEP::RandFlat   randFlat_;
      CLHEP::RandGaussQ randGaussQ_;
      std::string inputRootFile_;
      std::string treeName_;
      VDResampler::TrainState state_;
      // Accumulates pz-fallback hits over the TTree loop; one summary warning in endJob
      // (V2 slope basis only; see PzFallbackStats).
      VDResampler::PzFallbackStats pzFallback_;
  };

  VDResamplerTrainFromRoot::VDResamplerTrainFromRoot(const Parameters& conf) :
    art::EDAnalyzer(conf),
    engine_     (createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    randFlat_   (engine_),
    randGaussQ_ (engine_),
    inputRootFile_(conf().InputRootFile()),
    treeName_     (conf().TreeName())
  {
    // Populate shared state from fhicl config
    state_.allAtOnceModelFile  = conf().SBDMallAtOnceModelFile();
    state_.stage1ModelFile     = conf().SBDMstage1ModelFile();
    state_.stage2ModelFile     = conf().SBDMstage2ModelFile();
    state_.ckptAllAtOnceFile   = conf().SBDMloadCheckPointAllAtOnceModelFile();
    state_.ckptStage1File      = conf().SBDMloadCheckPointStage1ModelFile();
    state_.ckptStage2File      = conf().SBDMloadCheckPointStage2ModelFile();
    state_.useTwoStageTraining = conf().SBDMuseTwoStageTraining();
    state_.momentumBasis       = VDResampler::parseMomentumBasis(conf().SBDMmomentumBasis(), "VDResamplerTrainFromRoot");
    state_.virtualDetectorID   = conf().VirtualDetectorID();
    state_.VDz0                = conf().VDz0();
    state_.VDr                 = conf().VDr();
    state_.pdgID               = conf().pdgID();
    state_.trainingEpochs      = conf().SBDMtrainingEpochs();
    state_.trainingSize        = conf().SBDMtrainingSize();
    state_.saveEpochs          = conf().SaveEpochs();
    state_.autoPlanner              = conf().SBDMautoCurriculumPlanner();
    state_.plannerSmoothWindow      = conf().SBDMplannerSmoothWindow();
    state_.plannerMinDelta          = conf().SBDMplannerMinDelta();
    state_.plannerPatience          = conf().SBDMplannerPatience();
    state_.plannerMinEpochsPerPhase = conf().SBDMplannerMinEpochsPerPhase();
    state_.plannerMaxEpochsPerPhase = conf().SBDMplannerMaxEpochsPerPhase();
    state_.curriculumEpochs              = conf().SBDMtrainingCurriculumEpochs();
    state_.curriculumLossWeightPower     = conf().SBDMtrainingCurriculumLossWeightPower();
    state_.curriculumGradientClip        = conf().SBDMtrainingCurriculumGradientClip();
    state_.curriculumLearningRate        = conf().SBDMtrainingCurriculumLearningRate();
    state_.curriculumBiasLowSigma        = conf().SBDMtrainingCurriculumBiasLowSigma();
    state_.curriculumTLowBound           = conf().SBDMtrainingCurriculumTLowBound();
    state_.curriculumBatchSize           = conf().SBDMtrainingCurriculumBatchSize();
    state_.promoteEMAOnStart             = conf().SBDMpromoteEMA();
    state_.curriculumPromoteEMA          = conf().SBDMtrainingCurriculumPromoteEMA();
    state_.curriculumUseDimWeightController = conf().SBDMtrainingCurriculumUseDimWeightController();
    state_.curriculumUsePeakWindowLoss   = conf().SBDMtrainingCurriculumUsePeakWindowLoss();
    state_.curriculumUsePeakSampling     = conf().SBDMtrainingCurriculumUsePeakSampling();
    state_.curriculumPeakAlpha           = conf().SBDMtrainingCurriculumPeakAlpha();
    state_.curriculumTFocusLow           = conf().SBDMtrainingCurriculumTFocusLow();
    state_.curriculumTFocusHigh          = conf().SBDMtrainingCurriculumTFocusHigh();
    state_.curriculumTFocusFraction      = conf().SBDMtrainingCurriculumTFocusFraction();
    // Prefer the new SamplesDrawnPerEpoch key; fall back to the deprecated SubsetSizePerEpoch
    // name (warn once) so existing plan files keep working.
    if (!conf().SBDMtrainingCurriculumSamplesDrawnPerEpoch().empty()) {
        state_.curriculumSamplesDrawnPerEpoch = conf().SBDMtrainingCurriculumSamplesDrawnPerEpoch();
    } else {
        state_.curriculumSamplesDrawnPerEpoch = conf().SBDMtrainingCurriculumSubsetSizePerEpoch();
        if (!state_.curriculumSamplesDrawnPerEpoch.empty())
            mf::LogWarning("VDResamplerTrainFromRoot")
                << "SBDMtrainingCurriculumSubsetSizePerEpoch is deprecated; "
                << "rename it to SBDMtrainingCurriculumSamplesDrawnPerEpoch.";
    }
    state_.curriculumMinDelta            = conf().SBDMtrainingCurriculumMinDelta();
    state_.denoiseDiagnosticTs           = conf().SBDMdenoiseDiagnosticTs();
    state_.denoiseDiagnosticSamples      = conf().SBDMdenoiseDiagnosticSamples();
    state_.denoiseDiagnosticUseEMA       = conf().SBDMdenoiseDiagnosticUseEMA();
    state_.partialReverseT0s             = conf().SBDMpartialReverseT0s();
    state_.partialReverseSamples         = conf().SBDMpartialReverseSamples();
    state_.partialReverseUseHeun         = conf().SBDMpartialReverseUseHeun();
    state_.partialReverseUseSDE          = conf().SBDMpartialReverseUseSDE();
    state_.partialReverseDiffusionSteps  = conf().SBDMpartialReverseDiffusionSteps();
    state_.partialReverseSdeToOdeSigmaThreshold = conf().SBDMpartialReverseSdeToOdeSigmaThreshold();
    state_.condLossDiagnostic            = conf().SBDMcondLossDiagnostic();
    state_.condLossDiagnosticSamples     = conf().SBDMcondLossDiagnosticSamples();
    state_.featureBlockDiagnostic        = conf().SBDMfeatureBlockDiagnostic();
    state_.featureBlockDiagnosticSamples = conf().SBDMfeatureBlockDiagnosticSamples();
    VDResampler::assemblePeakWindows(state_,
        conf().SBDMpeakWindowDims(), conf().SBDMpeakWindowLows(), conf().SBDMpeakWindowHighs(),
        conf().SBDMpeakGMaxes(), conf().SBDMpeakSigma0s(), conf().SBDMpeakAlphas(), "VDResamplerTrainFromRoot");

    VDResampler::validateGeometry(state_.VDr, state_.VDz0, "VDResamplerTrainFromRoot");
    // Resolve the prediction target up front so the curriculum builder can coerce
    // lossWeightPower to 0 under V (and so p.predictionTarget reuses the same result).
    ScoreBasedDiffusionModel::PredictionTarget predictionTarget;
    {
        bool epsValue = false;
        bool epsSet = conf().SBDMepsPrediction(epsValue); // OptionalAtom: true if present
        predictionTarget = VDResampler::resolvePredictionTarget(
            conf().SBDMpredictionTarget(), epsSet, epsValue, "VDResamplerTrainFromRoot");
    }

    // Prefer the new SBDMsamplesDrawnPerEpoch; fall back to the deprecated
    // SBDMtrainingSubsetSizePerEpoch (warn once) so existing configs keep working. This is the
    // non-curriculum default; the per-phase curriculum key is handled separately above.
    int samplesDrawnPerEpoch = conf().SBDMsamplesDrawnPerEpoch();
    if (samplesDrawnPerEpoch == 0 && conf().SBDMtrainingSubsetSizePerEpoch() != 0) {
        samplesDrawnPerEpoch = conf().SBDMtrainingSubsetSizePerEpoch();
        mf::LogWarning("VDResamplerTrainFromRoot")
            << "SBDMtrainingSubsetSizePerEpoch is deprecated; rename it to SBDMsamplesDrawnPerEpoch.";
    }

    VDResampler::validateAndBuildCurriculum(state_, "VDResamplerTrainFromRoot",
        conf().SBDMlossWeightPower(), conf().SBDMgradientClip(), conf().SBDMlearningRate(),
        conf().SBDMbiasLowSigma(),   conf().SBDMtLowBound(),    conf().SBDMbatchSize(),
        conf().SBDMtFocusLow(),      conf().SBDMtFocusHigh(),   conf().SBDMtFocusFraction(),
        /*defaultPromoteEMA=*/false, conf().SBDMuseDimWeightController(),
        samplesDrawnPerEpoch, predictionTarget);

    VDResampler::ModelBuildParams p;
    p.timeEmbeddingDim           = conf().SBDMtimeEmbeddingDim();
    p.inputEmbeddingDims         = conf().SBDMinputEmbeddingDims();
    p.conditionEmbeddingDims     = conf().SBDMconditionEmbeddingDims();
    p.hidden                     = conf().SBDMhidden();
    p.layers                     = conf().SBDMlayers();
    p.optimizer                  = VDResampler::parseOptimizer(conf().SBDMoptimizer(), "VDResamplerTrainFromRoot");
    p.adamBeta1                  = conf().SBDMadamBeta1();
    p.adamBeta2                  = conf().SBDMadamBeta2();
    p.adamEps                    = conf().SBDMadamEps();
    p.noiseSchedule              = VDResampler::parseNoiseSchedule(conf().SBDMnoiseSchedule(), "VDResamplerTrainFromRoot");
    p.betaMin                    = conf().SBDMbetaMin();
    p.betaMax                    = conf().SBDMbetaMax();
    p.cosineOffset               = conf().SBDMcosineOffset();
    p.logSigMin                  = conf().SBDMlogSigMin();
    p.logSigMax                  = conf().SBDMlogSigMax();
    p.predictionTarget           = predictionTarget; // resolved above (before the curriculum build)
    p.lossWeightPower            = conf().SBDMlossWeightPower();
    p.batchSize                  = conf().SBDMbatchSize();
    p.gradientClip               = conf().SBDMgradientClip();
    p.learningRate               = conf().SBDMlearningRate();
    p.useDimWeightController     = conf().SBDMuseDimWeightController();
    p.dimWeightControllerEMADecay= conf().SBDMdimWeightControllerEMADecay();
    p.useEMANetwork              = conf().SBDMuseEMANetwork();
    p.emaNetworkDecay            = conf().SBDMemaNetworkDecay();
    p.diffusionSteps             = conf().SBDMdiffusionSteps();
    VDResampler::buildModels(state_, p, randFlat_, randGaussQ_, "VDResamplerTrainFromRoot");
  }

  void VDResamplerTrainFromRoot::endJob() {
    // Open ROOT file and read training data. TFile::Open is used to handle xroot:// paths.
    auto fin = std::unique_ptr<TFile>{TFile::Open(inputRootFile_.c_str(), "READ")};
    if (!fin || fin->IsZombie())
      throw cet::exception("VDResamplerTrainFromRoot") << "Cannot open ROOT file: " << inputRootFile_;
    TTree* ttree = dynamic_cast<TTree*>(fin->Get(treeName_.c_str()));
    if (!ttree)
      throw cet::exception("VDResamplerTrainFromRoot") << "Cannot find TTree: " << treeName_;

    double time, x, y, z, px, py, pz;
    int stepPdgId;
    ULong64_t vdId;
    ttree->SetBranchAddress("time",            &time);
    ttree->SetBranchAddress("x",               &x);
    ttree->SetBranchAddress("y",               &y);
    ttree->SetBranchAddress("z",               &z);
    ttree->SetBranchAddress("px",              &px);
    ttree->SetBranchAddress("py",              &py);
    ttree->SetBranchAddress("pz",              &pz);
    ttree->SetBranchAddress("pdgId",           &stepPdgId);
    ttree->SetBranchAddress("virtualdetectorId", &vdId);

    for (Long64_t i = 0; i < ttree->GetEntries(); ++i) {
      ttree->GetEntry(i);
      if (vdId != state_.virtualDetectorID || (stepPdgId != state_.pdgID && state_.pdgID != 0) || pz <= 0)
        continue;

      // mom0/mom1/mom2 are the three transformed momentum values in state_.momentumBasis.
      double x_trans, y_trans, t_trans, mom0_t, mom1_t, mom2_t;
      VDResampler::forwardTransformSample(
          x, y, z, time, px, py, pz,
          VDResampler::kX0, VDResampler::kY0, VDResampler::kT0, VDResampler::kTScale, VDResampler::kP0,
          state_.VDr, state_.VDz0,
          x_trans, y_trans, t_trans, mom0_t, mom1_t, mom2_t,
          state_.momentumBasis, &pzFallback_);

      VDResampler::accumulateNorm(state_, t_trans, x_trans, y_trans, mom0_t, mom1_t, mom2_t);
      VDResampler::collectSample (state_, t_trans, x_trans, y_trans, mom0_t, mom1_t, mom2_t);
    }
    fin->Close();

    // Single summary warning for any pz fallbacks accumulated over the TTree loop.
    if (pzFallback_.count > 0) {
      std::ostringstream oss;
      oss << "pz fell below kPzSafetyEpsilon (" << VDResampler::kPzSafetyEpsilon << ") in "
          << pzFallback_.count << " hit(s); the floor was used in a slope division "
          << "(pz>0 expected from the selection). First " << pzFallback_.firstValues.size()
          << " offending pz value(s):";
      for (double v : pzFallback_.firstValues) oss << ' ' << v;
      mf::LogWarning("VDResamplerTrainFromRoot") << oss.str();
    }

    VDResampler::runTraining(state_, "VDResamplerTrainFromRoot");
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::VDResamplerTrainFromRoot)
