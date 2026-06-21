// VDResamplerTrain_module.cc
// For StepPointMCs on the designated virtual detectors, train either:
//   1) a two-stage score-based diffusion model with
//      Stage 1: unconditional model for (t', x', y')
//      Stage 2: conditional model for (p_r', p_phi', p_z' | t', x', y')
//   2) a single unconditional score-based diffusion model for
//      (t', x', y', p_r', p_phi', p_z')
// and store the trained model parameters in CSV files.
// note that p_z are filtered and only hits with positive p_z are kept
// Yongyi Wu, Mar. 2026

#include <memory>

#include "Offline/STMMC/inc/VDResamplerTrainCommon.hh"
#include "Offline/STMMC/inc/VDResamplerTransforms.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/SeedService/inc/SeedService.hh"

namespace mu2e {
  class VDResamplerTrain : public art::EDAnalyzer {
    public:
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> StepPointMCsTag{ Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs") };
        fhicl::Atom<art::InputTag> SimParticlemvTag{ Name("SimParticlemvTag"), Comment("Tag identifying the SimParticlemv") };
        fhicl::Atom<bool>        SBDMuseTwoStageTraining{    Name("SBDMuseTwoStageTraining"),    Comment("If true, train two-stage model; if false, train all 6 dimensions at once."), true };
        fhicl::Atom<std::string> SBDMallAtOnceModelFile{    Name("SBDMallAtOnceModelFile"),    Comment("CSV filename for the all-at-once 6D model"),  "" };
        fhicl::Atom<std::string> SBDMstage1ModelFile{        Name("SBDMstage1ModelFile"),        Comment("CSV filename for the stage-1 model"),         "" };
        fhicl::Atom<std::string> SBDMstage2ModelFile{        Name("SBDMstage2ModelFile"),        Comment("CSV filename for the stage-2 model"),         "" };
        fhicl::Atom<std::string> SBDMloadCheckPointAllAtOnceModelFile{ Name("SBDMloadCheckPointAllAtOnceModelFile"), Comment("Load checkpoint CSV for all-at-once model"), "" };
        fhicl::Atom<std::string> SBDMloadCheckPointStage1ModelFile{    Name("SBDMloadCheckPointStage1ModelFile"),    Comment("Load checkpoint CSV for stage-1 model"),    "" };
        fhicl::Atom<std::string> SBDMloadCheckPointStage2ModelFile{    Name("SBDMloadCheckPointStage2ModelFile"),    Comment("Load checkpoint CSV for stage-2 model"),    "" };
        fhicl::Atom<bool>       SBDMpromoteEMA{                           Name("SBDMpromoteEMA"),                           Comment("Promote EMA weights to network (and reset optimizer) once at the start of training"), false };
        fhicl::Atom<int>    VirtualDetectorID{ Name("VirtualDetectorID"), Comment("ID of the virtual detector to train on"),    116 };
        fhicl::Atom<double> VDz0{              Name("VDz0"),              Comment("z coordinate of the virtual detector"),      37700.39 };
        fhicl::Atom<double> VDr{               Name("VDr"),               Comment("VD radius"),                                 2000.0 };
        fhicl::Atom<int>    pdgID{             Name("pdgID"),             Comment("pdgID of the particle to train on"),         22 };
        fhicl::Atom<int>    SBDMtimeEmbeddingDim{ Name("SBDMtimeEmbeddingDim"), Comment("Time embedding dimension"),            0 };
        fhicl::Sequence<int> SBDMinputEmbeddingDims{     Name("SBDMinputEmbeddingDims"),     Comment("Per-state-dim Fourier depth: [] none, [k] broadcast, or length-dim list; each 0 or even >= 2"),     std::vector<int>() };
        fhicl::Sequence<int> SBDMconditionEmbeddingDims{ Name("SBDMconditionEmbeddingDims"), Comment("Per-condition-dim Fourier depth: [] none, [k] broadcast, or length-condDim list; each 0 or even >= 2"), std::vector<int>() };
        fhicl::Atom<int>    SBDMhidden{        Name("SBDMhidden"),        Comment("Size of hidden layers"),                     128 };
        fhicl::Atom<int>    SBDMlayers{        Name("SBDMlayers"),        Comment("Number of layers"),                          4 };
        fhicl::Atom<std::string> SBDMoptimizer{ Name("SBDMoptimizer"),   Comment("Optimizer (SGD or ADAM)"),                   "ADAM" };
        fhicl::Atom<double> SBDMadamBeta1{     Name("SBDMadamBeta1"),     Comment("Adam beta1"),                                0.9 };
        fhicl::Atom<double> SBDMadamBeta2{     Name("SBDMadamBeta2"),     Comment("Adam beta2"),                                0.999 };
        fhicl::Atom<double> SBDMadamEps{       Name("SBDMadamEps"),       Comment("Adam epsilon"),                              1e-8 };
        fhicl::Atom<std::string> SBDMnoiseSchedule{ Name("SBDMnoiseSchedule"), Comment("Noise schedule (LINEAR/COSINE/LOGSIG)"), "COSINE" };
        fhicl::Atom<double> SBDMbetaMin{       Name("SBDMbetaMin"),       Comment("Min beta (LINEAR schedule)"),                1e-4 };
        fhicl::Atom<double> SBDMbetaMax{       Name("SBDMbetaMax"),       Comment("Max beta (LINEAR schedule)"),                0.02 };
        fhicl::Atom<double> SBDMcosineOffset{  Name("SBDMcosineOffset"),  Comment("Offset (COSINE schedule)"),                  0.008 };
        fhicl::Atom<double> SBDMlogSigMin{     Name("SBDMlogSigMin"),     Comment("Min sigma (LOGSIG schedule)"),               1e-5 };
        fhicl::Atom<double> SBDMlogSigMax{     Name("SBDMlogSigMax"),     Comment("Max sigma (LOGSIG schedule)"),               1.0 };
        fhicl::Atom<bool>   SBDMepsPrediction{ Name("SBDMepsPrediction"), Comment("Predict eps (true) or score (false)"),       false };
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
        fhicl::Atom<int>    SBDMtrainingSubsetSizePerEpoch{ Name("SBDMtrainingSubsetSizePerEpoch"), Comment("Subset size per epoch (0 = full set)"), 0 };
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
        fhicl::Sequence<double> SBDMtrainingCurriculumPeakAlpha{ Name("SBDMtrainingCurriculumPeakAlpha"), Comment("Per-phase override of peak importance-sampling alpha for ALL windows (1=unbiased, <1=up-weight feature); negative = use base SBDMpeakAlphas. Empty = no override."), std::vector<double>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumTFocusLow{          Name("SBDMtrainingCurriculumTFocusLow"),          Comment("t focus window low edge per curriculum phase"),       std::vector<double>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumTFocusHigh{         Name("SBDMtrainingCurriculumTFocusHigh"),         Comment("t focus window high edge per curriculum phase"),      std::vector<double>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumTFocusFraction{     Name("SBDMtrainingCurriculumTFocusFraction"),     Comment("Realized in-window t focus fraction per curriculum phase (rest drawn from complement)"),       std::vector<double>() };
        fhicl::Sequence<int>    SBDMtrainingCurriculumSubsetSizePerEpoch{ Name("SBDMtrainingCurriculumSubsetSizePerEpoch"), Comment("Training subset size per epoch per curriculum phase (0 = full set)"), std::vector<int>() };
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
      explicit VDResamplerTrain(const Parameters& conf);
      void analyze(const art::Event& e) override;
      void endJob() override;

    private:
      art::RandomNumberGenerator::base_engine_t& engine_;
      CLHEP::RandFlat   randFlat_;
      CLHEP::RandGaussQ randGaussQ_;
      art::ProductToken<StepPointMCCollection> StepPointMCsToken_;
      art::ProductToken<SimParticleCollection> SimParticlemvToken_;
      VDResampler::TrainState state_;
  };

  VDResamplerTrain::VDResamplerTrain(const Parameters& conf) :
    art::EDAnalyzer(conf),
    engine_     (createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    randFlat_   (engine_),
    randGaussQ_ (engine_),
    StepPointMCsToken_(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken_(consumes<SimParticleCollection>(conf().SimParticlemvTag()))
  {
    // Populate shared state from fhicl config
    state_.allAtOnceModelFile  = conf().SBDMallAtOnceModelFile();
    state_.stage1ModelFile     = conf().SBDMstage1ModelFile();
    state_.stage2ModelFile     = conf().SBDMstage2ModelFile();
    state_.ckptAllAtOnceFile   = conf().SBDMloadCheckPointAllAtOnceModelFile();
    state_.ckptStage1File      = conf().SBDMloadCheckPointStage1ModelFile();
    state_.ckptStage2File      = conf().SBDMloadCheckPointStage2ModelFile();
    state_.useTwoStageTraining = conf().SBDMuseTwoStageTraining();
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
    state_.curriculumPeakAlpha           = conf().SBDMtrainingCurriculumPeakAlpha();
    state_.curriculumTFocusLow           = conf().SBDMtrainingCurriculumTFocusLow();
    state_.curriculumTFocusHigh          = conf().SBDMtrainingCurriculumTFocusHigh();
    state_.curriculumTFocusFraction      = conf().SBDMtrainingCurriculumTFocusFraction();
    state_.curriculumSubsetSizePerEpoch  = conf().SBDMtrainingCurriculumSubsetSizePerEpoch();
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
        conf().SBDMpeakGMaxes(), conf().SBDMpeakSigma0s(), conf().SBDMpeakAlphas(), "VDResamplerTrain");

    VDResampler::validateGeometry(state_.VDr, state_.VDz0, "VDResamplerTrain");
    VDResampler::validateAndBuildCurriculum(state_, "VDResamplerTrain",
        conf().SBDMlossWeightPower(), conf().SBDMgradientClip(), conf().SBDMlearningRate(),
        conf().SBDMbiasLowSigma(),   conf().SBDMtLowBound(),    conf().SBDMbatchSize(),
        conf().SBDMtFocusLow(),      conf().SBDMtFocusHigh(),   conf().SBDMtFocusFraction(),
        /*defaultPromoteEMA=*/false, conf().SBDMuseDimWeightController(),
        conf().SBDMtrainingSubsetSizePerEpoch());

    VDResampler::ModelBuildParams p;
    p.timeEmbeddingDim           = conf().SBDMtimeEmbeddingDim();
    p.inputEmbeddingDims         = conf().SBDMinputEmbeddingDims();
    p.conditionEmbeddingDims     = conf().SBDMconditionEmbeddingDims();
    p.hidden                     = conf().SBDMhidden();
    p.layers                     = conf().SBDMlayers();
    p.optimizer                  = VDResampler::parseOptimizer(conf().SBDMoptimizer(), "VDResamplerTrain");
    p.adamBeta1                  = conf().SBDMadamBeta1();
    p.adamBeta2                  = conf().SBDMadamBeta2();
    p.adamEps                    = conf().SBDMadamEps();
    p.noiseSchedule              = VDResampler::parseNoiseSchedule(conf().SBDMnoiseSchedule(), "VDResamplerTrain");
    p.betaMin                    = conf().SBDMbetaMin();
    p.betaMax                    = conf().SBDMbetaMax();
    p.cosineOffset               = conf().SBDMcosineOffset();
    p.logSigMin                  = conf().SBDMlogSigMin();
    p.logSigMax                  = conf().SBDMlogSigMax();
    p.epsPrediction              = conf().SBDMepsPrediction();
    p.lossWeightPower            = conf().SBDMlossWeightPower();
    p.batchSize                  = conf().SBDMbatchSize();
    p.gradientClip               = conf().SBDMgradientClip();
    p.learningRate               = conf().SBDMlearningRate();
    p.useDimWeightController     = conf().SBDMuseDimWeightController();
    p.dimWeightControllerEMADecay= conf().SBDMdimWeightControllerEMADecay();
    p.useEMANetwork              = conf().SBDMuseEMANetwork();
    p.emaNetworkDecay            = conf().SBDMemaNetworkDecay();
    p.diffusionSteps             = conf().SBDMdiffusionSteps();
    VDResampler::buildModels(state_, p, randFlat_, randGaussQ_, "VDResamplerTrain");
  }

  void VDResamplerTrain::analyze(const art::Event& event) {
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken_);
    if (StepPointMCs.empty()) return;
    auto const& SimParticles = event.getProduct(SimParticlemvToken_);
    if (SimParticles.empty()) return;

    for (const StepPointMC& step : StepPointMCs) {
      const SimParticle& particle = SimParticles.at(step.trackId());
      int    stepPdgId       = particle.pdgId();
      auto   vdId            = step.virtualDetectorId();
      double pz              = step.momentum().z();
      if (vdId != state_.virtualDetectorID || (stepPdgId != state_.pdgID && state_.pdgID != 0) || pz <= 0)
        continue;

      double x    = step.position().x();
      double y    = step.position().y();
      double z    = step.position().z();
      double px   = step.momentum().x();
      double py   = step.momentum().y();
      double time = step.time();

      double x_trans, y_trans, t_trans, pr_t, pphi_t, pz_t;
      VDResampler::forwardTransformSample(
          x, y, z, time, px, py, pz,
          VDResampler::kX0, VDResampler::kY0, VDResampler::kT0, VDResampler::kTScale, VDResampler::kP0,
          state_.VDr, state_.VDz0,
          x_trans, y_trans, t_trans, pr_t, pphi_t, pz_t);

      VDResampler::accumulateNorm(state_, t_trans, x_trans, y_trans, pr_t, pphi_t, pz_t);
      VDResampler::collectSample (state_, t_trans, x_trans, y_trans, pr_t, pphi_t, pz_t);
    }
  }

  void VDResamplerTrain::endJob() {
    VDResampler::runTraining(state_, "VDResamplerTrain");
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::VDResamplerTrain)
