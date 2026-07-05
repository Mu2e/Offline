#pragma once

// Shared state and logic for VDResamplerTrain and VDResamplerTrainFromRoot modules.
// Header-only — all functions are inline. Include pattern follows VDResamplerTransforms.hh.
// Yongyi Wu, Jun. 2026

#include <algorithm>
#include <cmath>
#include <deque>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "Offline/MachineLearningTools/inc/ScoreBasedDiffusionModel.hh"
#include "Offline/STMMC/inc/VDResamplerTransforms.hh"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"

namespace mu2e {
namespace VDResampler {

// ---------------------------------------------------------------------------
// TrainState — all mutable state shared between the two Train modules
// ---------------------------------------------------------------------------
struct TrainState {
    // Models
    std::unique_ptr<ScoreBasedDiffusionModel> allAtOnceModel;
    std::unique_ptr<ScoreBasedDiffusionModel> stage1Model;
    std::unique_ptr<ScoreBasedDiffusionModel> stage2Model;

    // Training data
    std::vector<DiffusionTrainingSample> allAtOnceTrainingData;
    std::vector<DiffusionTrainingSample> stage1TrainingData;
    std::vector<DiffusionTrainingSample> stage2TrainingData;

    // Output file paths
    std::string allAtOnceModelFile;
    std::string stage1ModelFile;
    std::string stage2ModelFile;

    // Checkpoint input file paths
    std::string ckptAllAtOnceFile;
    std::string ckptStage1File;
    std::string ckptStage2File;

    // Flags
    bool useTwoStageTraining  = true;
    bool trainStage1          = false;
    bool trainStage2          = false;

    // Geometry / particle filter
    unsigned long virtualDetectorID = 0;
    double VDz0 = 0.0;
    double VDr  = 0.0;
    int pdgID   = 0;

    // Training control
    int trainingEpochs              = 0;
    int trainingSize                = -1;
    std::vector<int> saveEpochs;
    bool saveAlsoCsv                = false; // if true, also write a CSV alongside every binary save

    // Automatic curriculum planner. When autoPlanner is true, each curriculum phase
    // trains until its (smoothed) loss converges instead of for a fixed epoch count;
    // see ConvergenceTracer and the planner branch of runTraining()'s trainLoop.
    bool   autoPlanner              = false;
    int    plannerSmoothWindow      = 10;    // trailing moving-average window W
    double plannerMinDelta          = 0.005; // relative improvement threshold (0.5%)
    int    plannerPatience          = 20;    // no-improvement epochs => converged
    int    plannerMinEpochsPerPhase = 30;    // floor; auto-raised to >= window+patience
    int    plannerMaxEpochsPerPhase = 100;   // hard cap; hitting it aborts with a warning

    // Curriculum
    int nPhase = 1;
    std::vector<int>    curriculumEpochs;
    std::vector<double> curriculumLossWeightPower;
    std::vector<double> curriculumGradientClip;
    std::vector<double> curriculumLearningRate;
    std::vector<bool>   curriculumBiasLowSigma;
    std::vector<double> curriculumTLowBound;
    std::vector<int>    curriculumBatchSize;
    std::vector<bool>   curriculumPromoteEMA;
    std::vector<bool>   curriculumUseDimWeightController;
    std::vector<bool>   curriculumUsePeakWindowLoss; // per-phase: plateau on peak-window loss (auto planner only)
    std::vector<bool>   curriculumUsePeakSampling;   // per-phase: enable peak emphasis sampling (empty => true all phases)
    std::vector<double> curriculumPeakAlpha;        // per-phase peak alpha override for ALL windows (<0 = use base)
    std::vector<double> curriculumTFocusLow;
    std::vector<double> curriculumTFocusHigh;
    std::vector<double> curriculumTFocusFraction;
    std::vector<int>    curriculumSubsetSizePerEpoch;
    std::vector<double> curriculumMinDelta; // per-phase planner minDelta (auto planner only)
    bool   promoteEMAOnStart   = false;
    bool   currentBiasLowSigma = false;
    double currentTLowBound    = 0.0;
    double currentTFocusLow      = 0.0;
    double currentTFocusHigh     = 0.0;
    double currentTFocusFraction = 0.0;
    int    currentSubsetSizePerEpoch = 0; // live trainingSubsetSizePerEpoch (per-phase under curriculum)
    bool   currentUsePeakWindowLoss = false; // live: track peak-window loss this phase (auto planner)
    std::vector<PeakWindow> currentPeakWindows; // live: s.peakWindows if phase enables sampling, else empty
    std::vector<int> phaseBoundaries;

    // One-step denoising diagnostic. When denoiseDiagnosticTs is non-empty, runTraining()
    // writes a denoising-diagnostic ROOT file for each active model INSTEAD of training it.
    std::vector<double> denoiseDiagnosticTs;
    int denoiseDiagnosticSamples = 100000; // samples per t value
    bool denoiseDiagnosticUseEMA = true;   // network choice for BOTH diagnostics; match the
                                           // generation config's useEMANetworkIfAvailable

    // Partial-reverse diagnostic. When partialReverseT0s is non-empty, runTraining() noises
    // data samples to each t0 and runs the full reverse sampler from there INSTEAD of
    // training. Much more expensive than the one-step diagnostic: one full reverse loop
    // (~t0*steps forward passes) per sample per t0.
    std::vector<double> partialReverseT0s;
    int    partialReverseSamples = 20000;       // samples per t0 value
    bool   partialReverseUseHeun = true;
    bool   partialReverseUseSDE  = true;
    int    partialReverseDiffusionSteps = -1;   // -1 = model's configured value
    double partialReverseSdeToOdeSigmaThreshold = -1.0;

    // Conditional-loss diagnostic (L_P(sigma) vs L_Q(sigma)). When true, runTraining() writes the
    // per-sample eps-loss-vs-sigma ROOT file (split by peak window) INSTEAD of training.
    bool   condLossDiagnostic        = false;
    int    condLossDiagnosticSamples = 200000;
    // Feature-block diagnostic (first-layer weight/grad L2 by input block + per-output-dim loss).
    bool   featureBlockDiagnostic        = false;
    int    featureBlockDiagnosticSamples = 50000; // draws used to accumulate the gradient

    // Peak importance sampling windows (see ScoreBasedDiffusionModel::train). Assembled by the
    // train module from parallel fcl sequences; held constant across the whole job (all curriculum
    // phases). Empty => disabled. Window edges are in TRANSFORMED (pre-z-score) units (e.g.
    // pz_t = log(pz/p0)); the model converts them to z-score via normalizeCoord. Dim is the state
    // coordinate index: all-at-once {t,x,y,pr,pphi,pz} -> pz=5; stage-2 state {pr,pphi,pz} -> pz=2.
    std::vector<PeakWindow> peakWindows;
    // Base (configured) alpha per window, snapshotted by assemblePeakWindows so a per-phase
    // override (curriculumPeakAlpha) can be applied non-destructively and the sentinel (<0)
    // can restore the original value. Parallel to peakWindows.
    std::vector<double> basePeakAlphas;

    // Welford normalization accumulators for the transformed coordinates. Slots
    // 0,1,2 are always (t, x, y); slots mom0,mom1,mom2 are the three momentum
    // values in the active MomentumBasis:
    //   V1_CylindricalTransformed : mom0=asinh(pr/p0), mom1=asinh(pphi/p0), mom2=log(pz/p0)
    //   V2_PtotSlopes             : mom0=log(pTotal/p0), mom1=ur,                 mom2=uphi
    //   V2_PtotSlopesAsinh        : mom0=log(pTotal/p0), mom1=asinh(ur/kUrSlopeScale), mom2=asinh(uphi/kUphiSlopeScale)
    double t_mean=0,    t_M2=0,    t_stdev=0;
    double x_mean=0,    x_M2=0,    x_stdev=0;
    double y_mean=0,    y_M2=0,    y_stdev=0;
    double mom0_mean=0, mom0_M2=0, mom0_stdev=0;
    double mom1_mean=0, mom1_M2=0, mom1_stdev=0;
    double mom2_mean=0, mom2_M2=0, mom2_stdev=0;
    int nNorm = 0;

    // Active momentum basis (set from fcl by the train modules; default V2).
    MomentumBasis momentumBasis = MomentumBasis::V2_PtotSlopes;

    // True for any V2 (pTotal/slopes) basis. V2 changes both the momentum slot
    // meaning and the two-stage split (stage1 = 1-D pTotal, stage2 = 5-D rest).
    bool isV2() const { return momentumBasis != MomentumBasis::V1_CylindricalTransformed; }
};

// ---------------------------------------------------------------------------
// ConvergenceTracer — detects "stable loss" within a curriculum phase.
//
// Strategy: patience on a smoothed loss. Per-epoch diffusion loss is noisy
// (stochastic t sampling + per-epoch subset sampling), so we smooth it with a
// trailing moving average over `window` epochs and watch the best smoothed value.
// A new smoothed value counts as an improvement only if it beats the best by a
// RELATIVE margin `minDelta`; after `patience` consecutive non-improving epochs
// (and at least `minEpochs` epochs into the phase) the phase is declared converged.
// reset() is called at every phase entry, so its state is independent of any epoch
// history a loaded checkpoint may carry in the model.
// ---------------------------------------------------------------------------
struct ConvergenceTracer {
    int    window;
    double minDelta;
    int    patience;
    int    minEpochs;

    std::deque<double> recent;  // last `window` raw losses
    double bestSmoothed = std::numeric_limits<double>::infinity();
    double bestRawLoss  = std::numeric_limits<double>::infinity();
    int    bestRawEpoch      = -1; // base-1 epoch index (within phase) of bestRawLoss
    int    bestSmoothedEpoch = -1; // base-1 epoch index (within phase) of bestSmoothed
    int    sinceImprove = 0;
    int    epochsThisPhase = 0;
    // Flags describing the most recent update() call; the planner uses them to decide
    // when to write the raw/smoothed best checkpoints and snapshot the network.
    bool   rawImproved      = false;
    bool   smoothedImproved = false;
    // Diagnostics from the most recent update(), for per-epoch planner logging.
    double lastSmoothed = std::numeric_limits<double>::quiet_NaN(); // current smoothed loss
    double lastDelta    = std::numeric_limits<double>::quiet_NaN(); // relative improvement vs prior best (>0 = better)

    void reset() {
        recent.clear();
        bestSmoothed = std::numeric_limits<double>::infinity();
        bestRawLoss  = std::numeric_limits<double>::infinity();
        bestRawEpoch      = -1;
        bestSmoothedEpoch = -1;
        sinceImprove = 0;
        epochsThisPhase = 0;
        rawImproved      = false;
        smoothedImproved = false;
        lastSmoothed = std::numeric_limits<double>::quiet_NaN();
        lastDelta    = std::numeric_limits<double>::quiet_NaN();
    }

    // Push one epoch's loss; returns true when the phase is considered converged.
    bool update(double loss) {
        ++epochsThisPhase;
        rawImproved      = false;
        smoothedImproved = false;
        if (loss < bestRawLoss) { bestRawLoss = loss; bestRawEpoch = epochsThisPhase; rawImproved = true; }

        recent.push_back(loss);
        if (static_cast<int>(recent.size()) > window) recent.pop_front();
        double smoothed = std::accumulate(recent.begin(), recent.end(), 0.0)
                          / static_cast<double>(recent.size());

        // Relative improvement of the smoothed loss vs the best seen so far (before this
        // update). Positive = improvement; NaN on the first epoch (no prior best yet).
        const double prevBest = bestSmoothed;
        lastSmoothed = smoothed;
        lastDelta    = std::isfinite(prevBest) ? (prevBest - smoothed) / prevBest
                                               : std::numeric_limits<double>::quiet_NaN();

        if (smoothed < bestSmoothed * (1.0 - minDelta)) {
            bestSmoothed = smoothed;
            bestSmoothedEpoch = epochsThisPhase;
            sinceImprove = 0;
            smoothedImproved = true;
        } else {
            ++sinceImprove;
        }
        return epochsThisPhase >= minEpochs && sinceImprove >= patience;
    }
};

// ---------------------------------------------------------------------------
// ModelBuildParams — SBDM constructor arguments, filled from fhicl Config
// ---------------------------------------------------------------------------
struct ModelBuildParams {
    int    timeEmbeddingDim = 0;
    std::vector<int> inputEmbeddingDims;     // per-state-dim Fourier depth ({} / {k} / length-dim)
    std::vector<int> conditionEmbeddingDims; // per-condition-dim Fourier depth ({} / {k} / length-condDim)
    int    hidden           = 128;
    int    layers           = 4;
    ScoreBasedDiffusionModel::OptimizerType     optimizer;
    double adamBeta1        = 0.9;
    double adamBeta2        = 0.999;
    double adamEps          = 1e-8;
    ScoreBasedDiffusionModel::NoiseScheduleType noiseSchedule;
    double betaMin          = 1e-4;
    double betaMax          = 0.02;
    double cosineOffset     = 0.008;
    double logSigMin        = 1e-5;
    double logSigMax        = 1.0;
    ScoreBasedDiffusionModel::PredictionTarget predictionTarget =
        ScoreBasedDiffusionModel::PredictionTarget::SCORE;
    double lossWeightPower  = 2.0;
    int    batchSize        = 32;
    double gradientClip     = 1.0;
    double learningRate     = 1e-3;
    bool   useDimWeightController        = false;
    double dimWeightControllerEMADecay   = 0.99;
    bool   useEMANetwork                 = false;
    double emaNetworkDecay               = 0.9999;
    int    diffusionSteps                = 200;
};

// ---------------------------------------------------------------------------
// parseOptimizer — convert fhicl string to enum
// ---------------------------------------------------------------------------
inline ScoreBasedDiffusionModel::OptimizerType
parseOptimizer(const std::string& opt, const std::string& moduleName) {
    if (opt == "SGD") return ScoreBasedDiffusionModel::OptimizerType::SGD;
    if (opt != "ADAM")
        mf::LogWarning(moduleName) << "Unrecognized SBDMoptimizer value \"" << opt << "\"; falling back to ADAM.";
    return ScoreBasedDiffusionModel::OptimizerType::ADAM;
}

// ---------------------------------------------------------------------------
// parseNoiseSchedule — convert fhicl string to enum
// ---------------------------------------------------------------------------
inline ScoreBasedDiffusionModel::NoiseScheduleType
parseNoiseSchedule(const std::string& sched, const std::string& moduleName) {
    if (sched == "LINEAR") return ScoreBasedDiffusionModel::NoiseScheduleType::LINEAR;
    if (sched == "LOGSIG") return ScoreBasedDiffusionModel::NoiseScheduleType::LOGSIG;
    if (sched != "COSINE")
        mf::LogWarning(moduleName) << "Unrecognized SBDMnoiseSchedule value \"" << sched << "\"; falling back to COSINE.";
    return ScoreBasedDiffusionModel::NoiseScheduleType::COSINE;
}

// ---------------------------------------------------------------------------
// parseMomentumBasis — convert fhicl string to MomentumBasis enum (if-chain, to
//   match parseOptimizer/parseNoiseSchedule; switch is not valid on std::string).
//   V1_CYLINDRICAL       : V1_CylindricalTransformed (legacy asinh(pr),asinh(pphi),log(pz))
//   V2_PTOT_SLOPES       : V2_PtotSlopes (log(pTotal), ur=pr/pz, uphi=pphi/pz) [default]
//   V2_PTOT_SLOPES_ASINH : V2_PtotSlopesAsinh (slopes wrapped in asinh)
// ---------------------------------------------------------------------------
inline MomentumBasis parseMomentumBasis(const std::string& b, const std::string& moduleName) {
    if (b == "V1_CYLINDRICAL")       return MomentumBasis::V1_CylindricalTransformed;
    if (b == "V2_PTOT_SLOPES")       return MomentumBasis::V2_PtotSlopes;
    if (b == "V2_PTOT_SLOPES_ASINH") return MomentumBasis::V2_PtotSlopesAsinh;
    mf::LogWarning(moduleName) << "Unrecognized SBDMmomentumBasis value \"" << b
                               << "\"; falling back to V2_PTOT_SLOPES.";
    return MomentumBasis::V2_PtotSlopes;
}

// ---------------------------------------------------------------------------
// ensureBinExtension — warn and fix a model file path that doesn't end in .bin
// ---------------------------------------------------------------------------
inline std::string ensureBinExtension(const std::string& path, const std::string& fieldName,
                                       const std::string& moduleName)
{
    if (path.empty()) return path;
    if (path.size() >= 4 && path.substr(path.size() - 4) == ".bin") return path;
    std::string fixed = path;
    if (fixed.size() >= 4 && fixed.substr(fixed.size() - 4) == ".csv")
        fixed = fixed.substr(0, fixed.size() - 4) + ".bin";
    else
        fixed += ".bin";
    mf::LogWarning(moduleName)
        << fieldName << " does not end with \".bin\" (got \"" << path
        << "\"); correcting to \"" << fixed << "\".";
    return fixed;
}

// ---------------------------------------------------------------------------
// validateGeometry — throw if VDr or VDz0 are geometrically invalid
// ---------------------------------------------------------------------------
inline void validateGeometry(double VDr, double VDz0, const std::string& moduleName) {
    if (VDr <= 0.0)
        throw cet::exception(moduleName)
            << "VDr must be positive (got " << VDr << "); "
            << "rho = r/VDr would produce inf/NaN in training data.";
    if (!std::isfinite(VDz0))
        throw cet::exception(moduleName) << "VDz0 must be finite (got " << VDz0 << ").";
}

// ---------------------------------------------------------------------------
// resolvePredictionTarget — pick the model's PredictionTarget from the new
//   SBDMpredictionTarget string and the deprecated optional SBDMepsPrediction bool.
//   Precedence: if the deprecated bool is set, warn; if BOTH are set, SBDMpredictionTarget
//   wins (warn on conflict). Otherwise the string (default "SCORE") is used.
// ---------------------------------------------------------------------------
inline ScoreBasedDiffusionModel::PredictionTarget resolvePredictionTarget(
    const std::string& targetStr, bool epsSet, bool epsValue, const std::string& moduleName)
{
    auto fromStr = ScoreBasedDiffusionModel::predictionTargetFromName(targetStr); // throws on bad string
    if (epsSet) {
        auto fromEps = epsValue ? ScoreBasedDiffusionModel::PredictionTarget::EPS
                                : ScoreBasedDiffusionModel::PredictionTarget::SCORE;
        if (fromStr != fromEps)
            mf::LogWarning(moduleName)
                << "Both SBDMpredictionTarget (" << targetStr << ") and the deprecated "
                << "SBDMepsPrediction (" << (epsValue ? "true" : "false")
                << ") are set and disagree; using SBDMpredictionTarget.";
        else
            mf::LogWarning(moduleName)
                << "SBDMepsPrediction is deprecated; use SBDMpredictionTarget instead.";
    }
    return fromStr;
}

// ---------------------------------------------------------------------------
// assemblePeakWindows — build s.peakWindows from the six parallel fcl sequences
//   (SBDMpeakWindowDims / Lows / Highs / GMaxes / Sigma0s / Alphas). All must have the same
//   length K (the number of windows); empty => peak sampling disabled. Per-element value
//   validation (dim range, gMax in (0,1), sum gMax < 1, ...) is performed by
//   ScoreBasedDiffusionModel::train(); here we only enforce equal lengths and copy into structs.
// ---------------------------------------------------------------------------
inline void assemblePeakWindows(TrainState& s,
    const std::vector<int>& dims, const std::vector<double>& lows, const std::vector<double>& highs,
    const std::vector<double>& gmaxes, const std::vector<double>& sigma0s, const std::vector<double>& alphas,
    const std::string& moduleName)
{
    const size_t K = dims.size();
    auto checkLen = [&](size_t n, const char* name) {
        if (n != K)
            throw cet::exception(moduleName)
                << "Peak window sequence length mismatch: " << name << " has " << n << ", expected "
                << K << " (all SBDMpeak* sequences must match SBDMpeakWindowDims).";
    };
    checkLen(lows.size(),    "SBDMpeakWindowLows");
    checkLen(highs.size(),   "SBDMpeakWindowHighs");
    checkLen(gmaxes.size(),  "SBDMpeakGMaxes");
    checkLen(sigma0s.size(), "SBDMpeakSigma0s");
    checkLen(alphas.size(),  "SBDMpeakAlphas");

    s.peakWindows.clear();
    s.basePeakAlphas.clear();
    for (size_t k = 0; k < K; ++k) {
        PeakWindow pw;
        pw.dim = dims[k]; pw.low = lows[k]; pw.high = highs[k];
        pw.gMax = gmaxes[k]; pw.sigma0 = sigma0s[k]; pw.alpha = alphas[k];
        s.peakWindows.push_back(pw);
        s.basePeakAlphas.push_back(alphas[k]); // snapshot for per-phase override restore
    }
    if (K > 0)
        mf::LogInfo(moduleName) << "Configured " << K << " peak importance-sampling window(s).";
}

// ---------------------------------------------------------------------------
// validateAndBuildCurriculum
//   Validates curriculum vectors in state, fills empty ones with defaults,
//   computes trainingEpochs, phaseBoundaries, and initialises currentBias/tLow.
//   No-op (except setting current phase-0 values) when curriculumEpochs is empty.
// ---------------------------------------------------------------------------
inline void validateAndBuildCurriculum(TrainState& s, const std::string& moduleName,
    double defaultLossWeightPower, double defaultGradientClip, double defaultLearningRate,
    bool defaultBiasLowSigma, double defaultTLowBound, int defaultBatchSize,
    double defaultTFocusLow = 0.0, double defaultTFocusHigh = 0.0, double defaultTFocusFraction = 0.0,
    bool defaultPromoteEMA = false, bool defaultUseDimWeightController = false,
    int defaultSubsetSizePerEpoch = 0,
    ScoreBasedDiffusionModel::PredictionTarget predictionTarget =
        ScoreBasedDiffusionModel::PredictionTarget::SCORE)
{
    // v-prediction's target already embeds the SNR (sigma) weighting, so lossWeightPower is
    // forced to 0 everywhere under V. Coerce BOTH the scalar broadcast default AND any
    // explicit per-phase array up front, so the schema table, the fill broadcast, and the
    // realized lossWeightPower_ all agree on 0.0 (updateLossWeightPower is the runtime backstop).
    if (predictionTarget == ScoreBasedDiffusionModel::PredictionTarget::V) {
        bool warned = (defaultLossWeightPower != 0.0);
        defaultLossWeightPower = 0.0;
        for (double& v : s.curriculumLossWeightPower) { if (v != 0.0) warned = true; v = 0.0; }
        if (warned)
            mf::LogWarning(moduleName)
                << "v-prediction selected: forcing all lossWeightPower (scalar and per-phase "
                << "curriculum) to 0.0 (the v target already embeds the SNR weighting).";
    }
    if (!s.curriculumEpochs.empty()) {
        s.nPhase = static_cast<int>(s.curriculumEpochs.size());
        if (s.nPhase == 1) {
            mf::LogWarning(moduleName) << "Only one curriculum phase. Curriculum training inputs will be ignored.";
        } else {
            s.trainingEpochs = std::accumulate(s.curriculumEpochs.begin(), s.curriculumEpochs.end(), 0);

            auto fill = [&](auto& vec, auto defVal, const char* name) {
                if (vec.empty()) {
                    vec.assign(s.nPhase, defVal);
                } else if (static_cast<int>(vec.size()) != s.nPhase) {
                    throw cet::exception(moduleName)
                        << "Inconsistent sizes for curriculum training parameters: " << name;
                }
            };
            fill(s.curriculumLossWeightPower, defaultLossWeightPower, "SBDMtrainingCurriculumLossWeightPower");
            fill(s.curriculumGradientClip,    defaultGradientClip,    "SBDMtrainingCurriculumGradientClip");
            fill(s.curriculumLearningRate,    defaultLearningRate,    "SBDMtrainingCurriculumLearningRate");
            fill(s.curriculumBiasLowSigma,    defaultBiasLowSigma,    "SBDMtrainingCurriculumBiasLowSigma");
            fill(s.curriculumTLowBound,       defaultTLowBound,       "SBDMtrainingCurriculumTLowBound");
            fill(s.curriculumBatchSize,       defaultBatchSize,       "SBDMtrainingCurriculumBatchSize");
            fill(s.curriculumPromoteEMA,      defaultPromoteEMA,      "SBDMtrainingCurriculumPromoteEMA");
            fill(s.curriculumUseDimWeightController, defaultUseDimWeightController, "SBDMtrainingCurriculumUseDimWeightController");
            fill(s.curriculumUsePeakWindowLoss, false, "SBDMtrainingCurriculumUsePeakWindowLoss");
            fill(s.curriculumUsePeakSampling, true, "SBDMtrainingCurriculumUsePeakSampling"); // default ON
            fill(s.curriculumPeakAlpha, -1.0, "SBDMtrainingCurriculumPeakAlpha"); // <0 sentinel => use base alpha
            fill(s.curriculumTFocusLow,       defaultTFocusLow,       "SBDMtrainingCurriculumTFocusLow");
            fill(s.curriculumTFocusHigh,      defaultTFocusHigh,      "SBDMtrainingCurriculumTFocusHigh");
            fill(s.curriculumTFocusFraction,  defaultTFocusFraction,  "SBDMtrainingCurriculumTFocusFraction");
            fill(s.curriculumSubsetSizePerEpoch, defaultSubsetSizePerEpoch, "SBDMtrainingCurriculumSubsetSizePerEpoch");
            fill(s.curriculumMinDelta,        s.plannerMinDelta,      "SBDMtrainingCurriculumMinDelta");

            // Validate focus-window values for every phase up front, so a bad late-phase
            // value fails at job start rather than at the phase boundary mid-run.
            for (int i = 0; i < s.nPhase; ++i) {
                double frac = s.curriculumTFocusFraction[i];
                if (frac < 0.0 || frac > 1.0)
                    throw cet::exception(moduleName)
                        << "SBDMtrainingCurriculumTFocusFraction[" << i << "] = " << frac
                        << " must be in [0,1]";
                if (frac > 0.0 &&
                    (s.curriculumTFocusLow[i] < 0.0 || s.curriculumTFocusHigh[i] > 1.0 ||
                     s.curriculumTFocusLow[i] >= s.curriculumTFocusHigh[i]))
                    throw cet::exception(moduleName)
                        << "Invalid t focus window in curriculum phase " << i
                        << ": require 0 <= low < high <= 1 (got ["
                        << s.curriculumTFocusLow[i] << ", " << s.curriculumTFocusHigh[i] << "])";
                // Per-phase peak-alpha override: negative = sentinel (use base); any other
                // value must be a valid alpha in [0,1] (matches train()'s per-window check).
                const double a = s.curriculumPeakAlpha[i];
                if (a >= 0.0 && a > 1.0)
                    throw cet::exception(moduleName)
                        << "SBDMtrainingCurriculumPeakAlpha[" << i << "] = " << a
                        << " must be in [0,1] (or negative to use the base SBDMpeakAlphas)";
            }
            // Warn if any phase requests a peak-alpha override but no windows are configured
            // (the override would be a no-op). Peak windows are global, so one check suffices.
            if (s.peakWindows.empty() &&
                std::any_of(s.curriculumPeakAlpha.begin(), s.curriculumPeakAlpha.end(),
                            [](double v) { return v >= 0.0; }))
                mf::LogWarning(moduleName)
                    << "SBDMtrainingCurriculumPeakAlpha sets an override for one or more phases "
                    << "but no peak windows are configured (SBDMpeakWindow* empty); the override "
                    << "has no effect.";
            // Warn if a phase disables peak sampling while no windows are configured anyway
            // (the toggle would be a no-op). Peak windows are global, so one check suffices.
            if (s.peakWindows.empty() &&
                std::any_of(s.curriculumUsePeakSampling.begin(), s.curriculumUsePeakSampling.end(),
                            [](bool v) { return !v; }))
                mf::LogWarning(moduleName)
                    << "SBDMtrainingCurriculumUsePeakSampling disables peak sampling for one or "
                    << "more phases but no peak windows are configured (SBDMpeakWindow* empty); "
                    << "the toggle has no effect.";

            s.phaseBoundaries.clear();
            int epochSum = 0;
            for (int e : s.curriculumEpochs) { epochSum += e; s.phaseBoundaries.push_back(epochSum); }

            std::stringstream ss;
            ss << "[Curriculum Training Schema]\n" << s.nPhase << " phases.\n";
            ss << std::setw(10) << "Epochs"
               << std::setw(20) << "Loss Weight Power"
               << std::setw(20) << "Gradient Clip"
               << std::setw(20) << "Learning Rate"
               << std::setw(15) << "BiasLowSigma"
               << std::setw(15) << "TLowBound"
               << std::setw(15) << "BatchSize"
               << std::setw(15) << "PromoteEMA"
               << std::setw(15) << "DimWeightCtl"
               << std::setw(15) << "TFocusLow"
               << std::setw(15) << "TFocusHigh"
               << std::setw(15) << "TFocusFrac"
               << std::setw(15) << "SubsetSize"
               << std::setw(15) << "MinDelta"
               << std::setw(15) << "PeakWinLoss"
               << std::setw(15) << "PeakAlpha" << "\n";
            for (int i = 0; i < s.nPhase; ++i)
                ss << std::setw(10) << s.curriculumEpochs[i]
                   << std::setw(20) << s.curriculumLossWeightPower[i]
                   << std::setw(20) << s.curriculumGradientClip[i]
                   << std::setw(20) << s.curriculumLearningRate[i]
                   << std::setw(15) << s.curriculumBiasLowSigma[i]
                   << std::setw(15) << s.curriculumTLowBound[i]
                   << std::setw(15) << s.curriculumBatchSize[i]
                   << std::setw(15) << s.curriculumPromoteEMA[i]
                   << std::setw(15) << s.curriculumUseDimWeightController[i]
                   << std::setw(15) << s.curriculumTFocusLow[i]
                   << std::setw(15) << s.curriculumTFocusHigh[i]
                   << std::setw(15) << s.curriculumTFocusFraction[i]
                   << std::setw(15) << s.curriculumSubsetSizePerEpoch[i]
                   << std::setw(15) << s.curriculumMinDelta[i]
                   << std::setw(15) << s.curriculumUsePeakWindowLoss[i]
                   << std::setw(15) << (s.curriculumPeakAlpha[i] >= 0.0
                                        ? std::to_string(s.curriculumPeakAlpha[i]) : std::string("base")) << "\n";
            ss << "[End of Curriculum Training Schema]";
            mf::LogInfo(moduleName) << ss.str();
        }
    }
    s.currentBiasLowSigma    = s.nPhase > 1 ? s.curriculumBiasLowSigma[0]    : defaultBiasLowSigma;
    s.currentTLowBound       = s.nPhase > 1 ? s.curriculumTLowBound[0]       : defaultTLowBound;
    s.currentTFocusLow       = s.nPhase > 1 ? s.curriculumTFocusLow[0]       : defaultTFocusLow;
    s.currentTFocusHigh      = s.nPhase > 1 ? s.curriculumTFocusHigh[0]      : defaultTFocusHigh;
    s.currentTFocusFraction  = s.nPhase > 1 ? s.curriculumTFocusFraction[0]  : defaultTFocusFraction;
    s.currentSubsetSizePerEpoch = s.nPhase > 1 ? s.curriculumSubsetSizePerEpoch[0] : defaultSubsetSizePerEpoch;
    s.currentUsePeakWindowLoss  = s.nPhase > 1 ? s.curriculumUsePeakWindowLoss[0]  : false;

    // Auto curriculum planner sanity checks. The min-epochs floor must cover both the
    // time for the moving average to fill (window) AND a full no-improvement streak
    // (patience) on top of it, otherwise convergence could fire spuriously early.
    if (s.autoPlanner) {
        if (s.plannerSmoothWindow < 1)
            throw cet::exception(moduleName)
                << "SBDMplannerSmoothWindow must be >= 1 (got " << s.plannerSmoothWindow << ")";
        if (s.plannerPatience < 1)
            throw cet::exception(moduleName)
                << "SBDMplannerPatience must be >= 1 (got " << s.plannerPatience << ")";
        if (s.plannerMaxEpochsPerPhase < 1)
            throw cet::exception(moduleName)
                << "SBDMplannerMaxEpochsPerPhase must be >= 1 (got " << s.plannerMaxEpochsPerPhase << ")";
        if (s.plannerMinDelta < 0.0 || s.plannerMinDelta >= 1.0)
            throw cet::exception(moduleName)
                << "SBDMplannerMinDelta must be in [0,1) (got " << s.plannerMinDelta << ")";
        if (s.nPhase > 1) {
            for (int i = 0; i < s.nPhase; ++i)
                if (s.curriculumMinDelta[i] < 0.0 || s.curriculumMinDelta[i] >= 1.0)
                    throw cet::exception(moduleName)
                        << "SBDMtrainingCurriculumMinDelta[" << i << "] = " << s.curriculumMinDelta[i]
                        << " must be in [0,1)";
        }
        int minFloor = s.plannerSmoothWindow + s.plannerPatience;
        if (s.plannerMinEpochsPerPhase < minFloor) {
            mf::LogWarning(moduleName)
                << "SBDMplannerMinEpochsPerPhase (" << s.plannerMinEpochsPerPhase
                << ") raised to " << minFloor << " to cover smoothWindow+patience.";
            s.plannerMinEpochsPerPhase = minFloor;
        }
        if (s.plannerMaxEpochsPerPhase < s.plannerMinEpochsPerPhase)
            throw cet::exception(moduleName)
                << "SBDMplannerMaxEpochsPerPhase (" << s.plannerMaxEpochsPerPhase
                << ") must be >= effective minEpochsPerPhase (" << s.plannerMinEpochsPerPhase << ")";

        // In planner mode each non-zero curriculum phase uses its configured epoch count
        // as that phase's max-epoch cap (overriding the global SBDMplannerMaxEpochsPerPhase).
        // A 0-epoch phase is a skipped "setup-only" phase. Any non-zero cap must leave
        // room to reach the convergence floor, else the phase could never converge and
        // would always abort.
        if (s.nPhase > 1) {
            for (int i = 0; i < s.nPhase; ++i) {
                if (s.curriculumEpochs[i] > 0 && s.curriculumEpochs[i] < s.plannerMinEpochsPerPhase)
                    throw cet::exception(moduleName)
                        << "Auto planner: SBDMtrainingCurriculumEpochs[" << i << "] = "
                        << s.curriculumEpochs[i] << " is used as that phase's max-epoch cap "
                        << "but is below the effective minEpochsPerPhase ("
                        << s.plannerMinEpochsPerPhase << "); the phase could never converge. "
                        << "Raise this phase's epochs to >= " << s.plannerMinEpochsPerPhase
                        << ", or lower SBDMplannerMinEpochsPerPhase / SBDMplannerSmoothWindow / "
                        << "SBDMplannerPatience.";
            }
        }

        // Peak-window plateau metric needs configured peak windows to produce a non-NaN
        // value. If any phase requests it but no windows exist, those phases silently fall
        // back to the aggregate loss every epoch — warn so the misconfig is visible at job
        // start. (Peak windows are global, so one check covers all phases.)
        if (s.peakWindows.empty() &&
            std::any_of(s.curriculumUsePeakWindowLoss.begin(), s.curriculumUsePeakWindowLoss.end(),
                        [](bool b) { return b; }))
            mf::LogWarning(moduleName)
                << "SBDMtrainingCurriculumUsePeakWindowLoss is true for one or more phases but no "
                << "peak windows are configured (SBDMpeakWindow* empty); those phases will fall "
                << "back to the aggregate loss (peak-window metric is NaN every epoch).";
    }
}

// Reject a checkpoint whose stored ModelLayout tag doesn't match the slot it is being
// loaded into, and ADOPT the checkpoint's momentum basis as authoritative for this run.
// The tag is packed into the SBDM's opaque basisTag() at save time (packBasisTag), so an
// all-at-once file dropped into a stage slot (or vice versa) loads self-consistently and
// would otherwise train/diagnose silently on the wrong architecture.
//
// Momentum basis handling: the training/diagnostic data is rebuilt from ROOT in
// s.momentumBasis, and the loaded network was trained in the checkpoint's basis. These MUST
// agree or the data (and the diagnostic true_dim values) live in a different momentum space
// than the network expects — e.g. raw ur/uphi vs asinh-wrapped slopes. Rather than force the
// user to keep SBDMmomentumBasis in sync with the checkpoint by hand, we take the basis FROM
// the checkpoint: s.momentumBasis is overwritten with the loaded tag's basis (warning if the
// configured value differed). buildModels runs before the ROOT event loop that consults
// s.momentumBasis, so this override governs the rebuilt data. When several checkpoints are
// loaded in one run they must agree; the second call detects a cross-stage disagreement.
//
// NOTE: a legacy v<=6 checkpoint predates the tag and reads back as basisTag()==0 ==
// (AllAtOnce6D, V1_CylindricalTransformed); feeding such a file into a stage slot is
// correctly flagged by the layout check here.
inline void checkModelLayout(const ScoreBasedDiffusionModel& model,
                             VDResampler::ModelLayout expected,
                             const std::string& file, const char* slot,
                             MomentumBasis& runBasis, bool& runBasisAdopted,
                             const std::string& moduleName) {
    const int tag = model.basisTag();
    mf::LogInfo(moduleName)
        << "Loaded " << slot << " checkpoint " << file << ": " << basisTagToString(tag);
    if (VDResampler::unpackModelLayout(tag) != expected)
        throw cet::exception(moduleName)
            << "Checkpoint " << file << " loaded into the " << slot
            << " slot has " << basisTagToString(tag) << ", but that slot requires layout "
            << modelLayoutName(expected)
            << ". The checkpoint parameter points at the wrong model file.";

    const MomentumBasis ckptBasis = VDResampler::unpackMomentumBasis(tag);
    if (!runBasisAdopted) {
        // First loaded checkpoint: adopt its basis for the whole run so the rebuilt data
        // matches the network. Warn if it differs from the configured SBDMmomentumBasis.
        if (ckptBasis != runBasis)
            mf::LogWarning(moduleName)
                << "Checkpoint " << file << " was trained in momentum basis "
                << momentumBasisName(ckptBasis) << ", which differs from the configured "
                << "SBDMmomentumBasis (" << momentumBasisName(runBasis)
                << "). Using the checkpoint's basis so the rebuilt data (and diagnostic "
                << "true_dim values) match the loaded network.";
        runBasis = ckptBasis;
        runBasisAdopted = true;
    } else if (ckptBasis != runBasis) {
        // A previously loaded checkpoint already fixed the run basis; this one disagrees.
        throw cet::exception(moduleName)
            << "Checkpoint " << file << " loaded into the " << slot << " slot has basis "
            << momentumBasisName(ckptBasis) << ", but another checkpoint in this run uses "
            << momentumBasisName(runBasis)
            << ". Stage checkpoints loaded together must share one momentum basis.";
    }
}

// ---------------------------------------------------------------------------
// buildModels — construct SBDM model objects (or load checkpoints) and reserve
//               training data vectors.  Sets trainStage1/trainStage2 flags.
// ---------------------------------------------------------------------------
inline void buildModels(TrainState& s, const ModelBuildParams& p,
                        CLHEP::RandFlat& rf, CLHEP::RandGaussQ& rg,
                        const std::string& moduleName)
{
    // Ensure all model output file paths end with .bin.  Checkpoint *load* paths
    // (ckpt*File) are left untouched: loadModel() auto-detects the format by
    // extension, so legacy .csv checkpoints must keep their original names.
    s.allAtOnceModelFile = ensureBinExtension(s.allAtOnceModelFile, "SBDMallAtOnceModelFile",  moduleName);
    s.stage1ModelFile    = ensureBinExtension(s.stage1ModelFile,    "SBDMstage1ModelFile",     moduleName);
    s.stage2ModelFile    = ensureBinExtension(s.stage2ModelFile,    "SBDMstage2ModelFile",     moduleName);

    // Mode/checkpoint consistency. Each training mode only consults the checkpoint
    // parameter(s) for the slot(s) it actually builds: all-at-once reads ckptAllAtOnceFile
    // and ignores the stage checkpoints; two-stage reads ckptStage1/2File and ignores the
    // all-at-once one. A checkpoint set on the wrong parameter is otherwise silently dropped
    // (no load, no error) and the model trains/diagnoses from random weights — an expensive,
    // invisible mistake. Reject the mismatch up front instead.
    if (s.useTwoStageTraining) {
        if (!s.ckptAllAtOnceFile.empty())
            throw cet::exception(moduleName)
                << "SBDMloadCheckPointAllAtOnceModelFile is set (" << s.ckptAllAtOnceFile
                << ") but SBDMuseTwoStageTraining=true, which never loads it. "
                << "Use SBDMloadCheckPointStage1ModelFile / SBDMloadCheckPointStage2ModelFile, "
                << "or set SBDMuseTwoStageTraining=false.";
    } else {
        if (!s.ckptStage1File.empty() || !s.ckptStage2File.empty())
            throw cet::exception(moduleName)
                << "A stage checkpoint is set (SBDMloadCheckPointStage1ModelFile=\""
                << s.ckptStage1File << "\", SBDMloadCheckPointStage2ModelFile=\""
                << s.ckptStage2File << "\") but SBDMuseTwoStageTraining=false, which never "
                << "loads them. Use SBDMloadCheckPointAllAtOnceModelFile, or set "
                << "SBDMuseTwoStageTraining=true.";
    }

    // Checkpoint-vs-slot-enabled consistency. A stage is active only when its output
    // SBDMstage{1,2}ModelFile is non-empty (that name gates trainStage{1,2} below and,
    // in diagnostic mode, becomes the output ROOT stem). If a checkpoint is set for a
    // stage whose model file is empty, the slot is skipped entirely — the checkpoint is
    // loaded-then-ignored (or never loaded), so the diagnostics/training the user asked
    // for silently never run. Reject that rather than no-op. (Two-stage only; the
    // all-at-once checkpoint is already handled by the mode check above.)
    if (s.useTwoStageTraining) {
        if (!s.ckptStage1File.empty() && s.stage1ModelFile.empty())
            throw cet::exception(moduleName)
                << "SBDMloadCheckPointStage1ModelFile is set (" << s.ckptStage1File
                << ") but SBDMstage1ModelFile is empty, so the stage-1 slot is disabled and "
                << "the checkpoint is never used. Set SBDMstage1ModelFile (it names the "
                << "output/diagnostic files), or clear the stage-1 checkpoint.";
        if (!s.ckptStage2File.empty() && s.stage2ModelFile.empty())
            throw cet::exception(moduleName)
                << "SBDMloadCheckPointStage2ModelFile is set (" << s.ckptStage2File
                << ") but SBDMstage2ModelFile is empty, so the stage-2 slot is disabled and "
                << "the checkpoint is never used. Set SBDMstage2ModelFile (it names the "
                << "output/diagnostic files), or clear the stage-2 checkpoint.";
    }

    // Validate EMA promotion requests against useEMANetwork
    bool anyPromotion = s.promoteEMAOnStart ||
        std::any_of(s.curriculumPromoteEMA.begin(), s.curriculumPromoteEMA.end(),
                    [](bool v){ return v; });
    if (anyPromotion && !p.useEMANetwork)
        throw cet::exception(moduleName)
            << "EMA promotion requested (SBDMpromoteEMA or SBDMtrainingCurriculumPromoteEMA) "
            << "but SBDMuseEMANetwork is false. Enable the EMA network to use this feature.";

    // Helper: build a fresh SBDM with phase-0 (or singleton) hyperparams
    auto makeSBDM = [&](int dim, int condDim) {
        double initLWP = s.nPhase > 1 ? s.curriculumLossWeightPower[0] : p.lossWeightPower;
        double initGC  = s.nPhase > 1 ? s.curriculumGradientClip[0]    : p.gradientClip;
        double initLR  = s.nPhase > 1 ? s.curriculumLearningRate[0]    : p.learningRate;
        bool   initUDWC= s.nPhase > 1 ? s.curriculumUseDimWeightController[0] : p.useDimWeightController;
        return std::make_unique<ScoreBasedDiffusionModel>(
            rf, rg, dim, condDim, p.timeEmbeddingDim, p.inputEmbeddingDims, p.conditionEmbeddingDims,
            p.hidden, p.layers, p.optimizer,
            p.adamBeta1, p.adamBeta2, p.adamEps, p.noiseSchedule,
            p.betaMin, p.betaMax, p.cosineOffset, p.logSigMin, p.logSigMax,
            p.predictionTarget, initLWP, p.batchSize, initGC, initLR,
            initUDWC, p.dimWeightControllerEMADecay,
            p.useEMANetwork, p.emaNetworkDecay, p.diffusionSteps
        );
    };

    // Phase-0 dim-weight-controller setting that should govern *this run*, regardless of
    // what a loaded checkpoint baked in. Applied as an explicit override after loadModel
    // (which otherwise restores the checkpoint's flag); updateUseDimWeightController logs
    // the frozen dimWeights when this turns the controller off.
    bool phase0UseDimWeightController =
        s.nPhase > 1 ? s.curriculumUseDimWeightController[0] : p.useDimWeightController;

    // Tracks whether a loaded checkpoint has fixed s.momentumBasis for this run yet.
    // checkModelLayout adopts the first loaded checkpoint's basis (overriding config) and
    // then enforces that any further loaded checkpoint agrees. Remains false when no
    // checkpoint is loaded (fresh training), leaving the configured basis in place.
    bool runBasisAdopted = false;

    if (s.useTwoStageTraining) {
        if (s.stage1ModelFile.empty() && s.stage2ModelFile.empty())
            throw cet::exception(moduleName)
                << "useTwoStageTraining requires at least one of SBDMstage1ModelFile and SBDMstage2ModelFile.";
        s.trainStage1 = !s.stage1ModelFile.empty();
        s.trainStage2 = !s.stage2ModelFile.empty();

        if (s.trainStage1) {
            if (!s.ckptStage1File.empty()) {
                s.stage1Model = std::make_unique<ScoreBasedDiffusionModel>(
                    ScoreBasedDiffusionModel::loadModel(rf, rg, s.ckptStage1File));
                checkModelLayout(*s.stage1Model,
                    s.isV2() ? VDResampler::ModelLayout::TwoStageStage1Ptot1D
                             : VDResampler::ModelLayout::AllAtOnce6D,
                    s.ckptStage1File, "stage-1", s.momentumBasis, runBasisAdopted, moduleName);
                s.stage1Model->updateUseDimWeightController(phase0UseDimWeightController);
                mf::LogInfo(moduleName)
                    << "Watch for parameter override: Loaded stage-1 model checkpoint from " << s.ckptStage1File;
            } else {
                // V1 two-stage: stage1 = 3-D (t,x,y). V2 two-stage: stage1 = 1-D (pTotal).
                s.stage1Model = makeSBDM(s.isV2() ? 1 : 3, 0);
            }
            if (s.trainingSize > 0) s.stage1TrainingData.reserve(s.trainingSize);
            else                    s.stage1TrainingData.reserve(1000);
        }
        if (s.trainStage2) {
            if (!s.ckptStage2File.empty()) {
                s.stage2Model = std::make_unique<ScoreBasedDiffusionModel>(
                    ScoreBasedDiffusionModel::loadModel(rf, rg, s.ckptStage2File));
                checkModelLayout(*s.stage2Model,
                    s.isV2() ? VDResampler::ModelLayout::TwoStageStage2_5D
                             : VDResampler::ModelLayout::AllAtOnce6D,
                    s.ckptStage2File, "stage-2", s.momentumBasis, runBasisAdopted, moduleName);
                s.stage2Model->updateUseDimWeightController(phase0UseDimWeightController);
                mf::LogInfo(moduleName)
                    << "Watch for parameter override: Loaded stage-2 model checkpoint from " << s.ckptStage2File;
            } else {
                // V1 stage2: 3-D (pr,pphi,pz) cond 3-D (t,x,y).
                // V2 stage2: 5-D (t,x,y,ur,uphi) cond 1-D (pTotal).
                s.stage2Model = s.isV2() ? makeSBDM(5, 1) : makeSBDM(3, 3);
            }
            if (s.trainingSize > 0) s.stage2TrainingData.reserve(s.trainingSize);
            else                    s.stage2TrainingData.reserve(1000);
        }
    } else {
        if (!s.ckptAllAtOnceFile.empty()) {
            s.allAtOnceModel = std::make_unique<ScoreBasedDiffusionModel>(
                ScoreBasedDiffusionModel::loadModel(rf, rg, s.ckptAllAtOnceFile));
            checkModelLayout(*s.allAtOnceModel, VDResampler::ModelLayout::AllAtOnce6D,
                s.ckptAllAtOnceFile, "all-at-once", s.momentumBasis, runBasisAdopted, moduleName);
            s.allAtOnceModel->updateUseDimWeightController(phase0UseDimWeightController);
            mf::LogInfo(moduleName)
                << "Watch for parameter override: Loaded all-at-once model checkpoint from " << s.ckptAllAtOnceFile;
        } else {
            s.allAtOnceModel = makeSBDM(6, 0);
        }
        if (s.trainingSize > 0) s.allAtOnceTrainingData.reserve(s.trainingSize);
        else                    s.allAtOnceTrainingData.reserve(1000);
    }
}

// ---------------------------------------------------------------------------
// accumulateNorm — Welford online update for normalization statistics
// ---------------------------------------------------------------------------
// mom0_t/mom1_t/mom2_t are the three transformed momentum values in s.momentumBasis
// (V1: asinh(pr/p0),asinh(pphi/p0),log(pz/p0); V2: log(pTotal/p0),ur,uphi).
inline void accumulateNorm(TrainState& s, double t_trans, double x_trans, double y_trans,
                           double mom0_t, double mom1_t, double mom2_t)
{
    if (s.nNorm >= s.trainingSize && s.trainingSize > 0) return;
    ++s.nNorm;
    auto update = [&](double& mean, double& M2, double val) {
        double delta = val - mean;
        mean += delta / static_cast<double>(s.nNorm);
        M2   += delta * (val - mean);
    };
    update(s.t_mean,    s.t_M2,    t_trans);
    update(s.x_mean,    s.x_M2,    x_trans);
    update(s.y_mean,    s.y_M2,    y_trans);
    update(s.mom0_mean, s.mom0_M2, mom0_t);
    update(s.mom1_mean, s.mom1_M2, mom1_t);
    update(s.mom2_mean, s.mom2_M2, mom2_t);
}

// ---------------------------------------------------------------------------
// collectSample — push a transformed hit into the appropriate training data vector(s).
//   The two-stage split is basis-dependent:
//     V1: stage1 = {t,x,y};               stage2 = {mom0,mom1,mom2} | {t,x,y}
//     V2: stage1 = {mom0=log(pTotal/p0)}; stage2 = {t,x,y,mom1=ur,mom2=uphi} | {mom0}
//   (full vector order: V1 (t,x,y,pr,pphi,pz); V2 (t,x,y,pTotal,ur,uphi).)
// ---------------------------------------------------------------------------
inline void collectSample(TrainState& s, double t_trans, double x_trans, double y_trans,
                           double mom0_t, double mom1_t, double mom2_t)
{
    if (s.useTwoStageTraining) {
        if (s.trainStage1) {
            DiffusionTrainingSample s1;
            // V2 stage1 is the 1-D pTotal model; V1 stage1 is the 3-D (t,x,y) model.
            if (s.isV2()) s1.x = {mom0_t};
            else          s1.x = {t_trans, x_trans, y_trans};
            s.stage1TrainingData.push_back(std::move(s1));
        }
        if (s.trainStage2) {
            DiffusionTrainingSample s2;
            if (s.isV2()) {
                s2.x    = {t_trans, x_trans, y_trans, mom1_t, mom2_t};
                s2.cond = {mom0_t};
            } else {
                s2.x    = {mom0_t, mom1_t, mom2_t};
                s2.cond = {t_trans, x_trans, y_trans};
            }
            s.stage2TrainingData.push_back(std::move(s2));
        }
    } else {
        DiffusionTrainingSample sa;
        sa.x = {t_trans, x_trans, y_trans, mom0_t, mom1_t, mom2_t};
        s.allAtOnceTrainingData.push_back(std::move(sa));
    }
}

// ---------------------------------------------------------------------------
// diagnosticTrueValueRanges — per-dimension [lo, hi] range (de-normalized, transformed units)
//   of the strided diagnostic sample set, padded, for histogram axes. Deterministic from the
//   data (independent of t/t0), so true and reconstructed values share one binning and overlay.
// ---------------------------------------------------------------------------
inline std::vector<std::pair<double,double>> diagnosticTrueValueRanges(
    const std::vector<DiffusionTrainingSample>& data, const std::vector<double>& mean,
    const std::vector<double>& stdev, int dim, size_t stride, size_t nUse)
{
    std::vector<std::pair<double,double>> r(dim, { std::numeric_limits<double>::max(),
                                                   std::numeric_limits<double>::lowest() });
    size_t used = 0;
    for (size_t k = 0; k < data.size() && used < nUse; k += stride, ++used)
        for (int i = 0; i < dim; ++i) {
            double v = data[k].x[i] * stdev[i] + mean[i];
            r[i].first  = std::min(r[i].first, v);
            r[i].second = std::max(r[i].second, v);
        }
    for (int i = 0; i < dim; ++i) {
        double lo = r[i].first, hi = r[i].second;
        if (!(hi > lo)) { lo -= 1.0; hi += 1.0; }       // guard against a degenerate (constant) dim
        double pad = 0.02 * (hi - lo);
        r[i] = { lo - pad, hi + pad };
    }
    return r;
}

// ---------------------------------------------------------------------------
// diagnosticBinCounts — per-dimension 1-D histogram bin count. Default `defaultBins` everywhere,
//   except dimensions targeted by a peak window: there the binning is made finer than the window's
//   sigma0 so the feature is resolved. sigma0 is in z-score units, the histogram axis in transformed
//   units, so the feature width on the axis is sigma0*stdev[dim]; we put kBinsPerSigma0 bins across
//   it. The smallest sigma0 wins if several windows share a dim. Capped to avoid pathological sizes.
// ---------------------------------------------------------------------------
inline std::vector<int> diagnosticBinCounts(
    int dim, int defaultBins, const std::vector<std::pair<double,double>>& ranges,
    const std::vector<double>& stdev, const std::vector<PeakWindow>& peakWindows)
{
    constexpr int kBinsPerSigma0 = 4;     // >= this many bins across one sigma0 (=> bin width < sigma0/4)
    constexpr int kMaxBins       = 20000;
    std::vector<int> nb(dim, defaultBins);
    for (const auto& pw : peakWindows) {
        if (pw.dim < 0 || pw.dim >= dim || pw.sigma0 <= 0.0) continue;
        const double featW = pw.sigma0 * stdev[pw.dim];               // feature width in transformed (axis) units
        const double span  = ranges[pw.dim].second - ranges[pw.dim].first;
        if (featW <= 0.0 || span <= 0.0) continue;
        const int needed = static_cast<int>(std::ceil(kBinsPerSigma0 * span / featW));
        nb[pw.dim] = std::min(kMaxBins, std::max(nb[pw.dim], needed));
    }
    return nb;
}

// ---------------------------------------------------------------------------
// runDenoiseDiagnostic — one-step denoising check of a (typically checkpoint-
//   loaded) model. For each requested diffusion time t, strides through the
//   normalized training data, noises each sample at t, runs a single forward
//   pass, and writes truth vs denoised values (in transformed, pre-z-score
//   units) to a ROOT file, one TTree per t value (named diag_t<value> with
//   '.' replaced by 'p', e.g. diag_t0p4). If a fine feature is present in the
//   denoised histogram but missing from generated samples, the reverse-process
//   sampler is the problem; if it is missing here too, the score network never
//   learned it. Branches: t, then per dimension i: true_dim<i>, denoised_dim<i>, noised_dim<i>
//   (the noised input the network sees), in the model's dimension order (all-at-once:
//   t,x,y,pr,pphi,pz; stage 2: pr,pphi,pz). Per t and dim it also writes histograms: truth/denoised/
//   noised 1-D overlays (*_h_true_dim<i>, *_h_denoised_dim<i>, *_h_noised_dim<i>; binning is finer on
//   any dim a peak window targets, resolving the window sigma0), a residual *_h_residual_dim<i>
//   (denoised-truth, auto-ranged), and a 2-D *_h2_true_vs_denoised_dim<i> (diagonal = perfect).
// ---------------------------------------------------------------------------
inline void runDenoiseDiagnostic(ScoreBasedDiffusionModel& model,
                                 const std::vector<DiffusionTrainingSample>& data, // already normalized
                                 const std::vector<double>& mean,
                                 const std::vector<double>& stdev,
                                 const TrainState& s,
                                 const std::string& outFile,
                                 const std::string& moduleName)
{
    TFile fout(outFile.c_str(), "RECREATE");
    if (fout.IsZombie())
        throw cet::exception(moduleName) << "Cannot open denoise diagnostic output file " << outFile;

    constexpr int kMaxDiagDims = 16;
    const int dim = static_cast<int>(data.front().x.size());
    if (dim > kMaxDiagDims)
        throw cet::exception(moduleName)
            << "Denoise diagnostic supports at most " << kMaxDiagDims << " dimensions (got " << dim << ")";

    const size_t nUse = (s.denoiseDiagnosticSamples > 0)
        ? std::min(data.size(), static_cast<size_t>(s.denoiseDiagnosticSamples))
        : data.size();
    const size_t stride = std::max<size_t>(1, data.size() / nUse);
    const auto ranges = diagnosticTrueValueRanges(data, mean, stdev, dim, stride, nUse);
    constexpr int kNValBins = 200;
    constexpr int kN2DCap   = 500; // 2-D histograms are at most this many bins per axis
    const auto nbins = diagnosticBinCounts(dim, kNValBins, ranges, stdev, s.peakWindows);
    // Auto-range the residual histograms (created with empty limits) via the buffer mechanism;
    // their scale varies with t. Restored before returning.
    const Int_t prevBufSize = TH1::GetDefaultBufferSize();
    TH1::SetDefaultBufferSize(50000);

    for (double t : s.denoiseDiagnosticTs) {
        std::ostringstream nm;
        nm << "diag_t" << t;
        std::string treeName = nm.str();
        std::replace(treeName.begin(), treeName.end(), '.', 'p');
        std::ostringstream title;
        title << "One-step denoising diagnostic at t=" << t;
        TTree* tree = new TTree(treeName.c_str(), title.str().c_str());

        double tBranch = t;
        double trueD[kMaxDiagDims]     = {0.0};
        double denoisedD[kMaxDiagDims] = {0.0};
        double noisedD[kMaxDiagDims]   = {0.0};
        tree->Branch("t", &tBranch);
        for (int i = 0; i < dim; ++i) {
            tree->Branch(("true_dim" + std::to_string(i)).c_str(), &trueD[i]);
            tree->Branch(("denoised_dim" + std::to_string(i)).c_str(), &denoisedD[i]);
            tree->Branch(("noised_dim" + std::to_string(i)).c_str(), &noisedD[i]);
        }
        // Per-dim histograms: truth / noised input / denoised overlay (shared binning, finer on peak
        // dims), residual (denoised-truth, auto-ranged), and 2-D truth-vs-denoised (diagonal = perfect).
        std::vector<TH1D*> hTrue(dim), hDen(dim), hNoised(dim), hResid(dim);
        std::vector<TH2D*> h2(dim);
        const std::string ts = std::to_string(t);
        for (int i = 0; i < dim; ++i) {
            const std::string di = std::to_string(i);
            const int n2d = std::min(kN2DCap, nbins[i]);
            hTrue[i] = new TH1D((treeName + "_h_true_dim" + di).c_str(),
                ("truth dim" + di + " at t=" + ts + ";value;count").c_str(),
                nbins[i], ranges[i].first, ranges[i].second);
            hDen[i] = new TH1D((treeName + "_h_denoised_dim" + di).c_str(),
                ("one-step denoised dim" + di + " at t=" + ts + ";value;count").c_str(),
                nbins[i], ranges[i].first, ranges[i].second);
            hNoised[i] = new TH1D((treeName + "_h_noised_dim" + di).c_str(),
                ("noised input dim" + di + " at t=" + ts + ";value;count").c_str(),
                nbins[i], ranges[i].first, ranges[i].second);
            hResid[i] = new TH1D((treeName + "_h_residual_dim" + di).c_str(),
                ("residual (denoised-truth) dim" + di + " at t=" + ts + ";denoised - truth;count").c_str(),
                nbins[i], 0.0, 0.0); // auto-ranged via buffer
            h2[i] = new TH2D((treeName + "_h2_true_vs_denoised_dim" + di).c_str(),
                ("truth vs denoised dim" + di + " at t=" + ts + ";truth;denoised").c_str(),
                n2d, ranges[i].first, ranges[i].second, n2d, ranges[i].first, ranges[i].second);
        }

        std::vector<double> noisedZ;
        size_t written = 0;
        for (size_t k = 0; k < data.size() && written < nUse; k += stride, ++written) {
            const auto& sample = data[k];
            auto res = model.denoiseOneStep(sample.x, sample.cond, t, s.denoiseDiagnosticUseEMA, &noisedZ);
            for (int i = 0; i < dim; ++i) {
                trueD[i]     = sample.x[i] * stdev[i] + mean[i];
                denoisedD[i] = res.value[i];
                noisedD[i]   = noisedZ[i] * stdev[i] + mean[i];
                hTrue[i]->Fill(trueD[i]);
                hDen[i]->Fill(denoisedD[i]);
                hNoised[i]->Fill(noisedD[i]);
                hResid[i]->Fill(denoisedD[i] - trueD[i]);
                h2[i]->Fill(trueD[i], denoisedD[i]);
            }
            tree->Fill();
        }
        tree->Write();
        for (int i = 0; i < dim; ++i) {
            hTrue[i]->Write(); hDen[i]->Write(); hNoised[i]->Write();
            hResid[i]->BufferEmpty(); // flush buffer so the auto-range is computed before writing
            hResid[i]->Write();
            h2[i]->Write();
        }
    }
    TH1::SetDefaultBufferSize(prevBufSize);
    fout.Close();
    mf::LogInfo(moduleName)
        << "Denoise diagnostic written to " << outFile
        << " (" << nUse << " samples per t, " << s.denoiseDiagnosticTs.size() << " t values)";
}

// ---------------------------------------------------------------------------
// runPartialReverseDiagnostic — multi-step reverse-sampler check of a (typically
//   checkpoint-loaded) model. For each requested start time t0, strides through
//   the normalized training data, noises each sample to t0, and runs the FULL
//   reverse sampler from t0 down to 0 (same integrator as generateSample).
//   Writes truth vs sampled values (transformed, pre-z-score units) to a ROOT
//   file, one TTree per t0 (named prev_t<value> with '.' replaced by 'p').
//   Scanning t0 localizes where generation loses a feature: if it survives a
//   start at t0 but not full generation, the t>t0 phase delivers a wrong
//   marginal; if it is already lost at t0, the score or its integration below
//   t0 is responsible. Branches: t0, then per dimension i: true_dim<i>, sampled_dim<i>,
//   noised_dim<i> (the noised t0 start state), in the model's dimension order. Per t0 and dim it
//   also writes histograms: truth/sampled/noised 1-D overlays (*_h_true_dim<i>, *_h_sampled_dim<i>,
//   *_h_noised_dim<i>; binning is finer on any dim a peak window targets), a residual
//   *_h_residual_dim<i> (sampled-truth, auto-ranged), and a 2-D *_h2_true_vs_sampled_dim<i>
//   (diagonal = perfect reproduction).
// ---------------------------------------------------------------------------
inline void runPartialReverseDiagnostic(ScoreBasedDiffusionModel& model,
                                        const std::vector<DiffusionTrainingSample>& data, // already normalized
                                        const std::vector<double>& mean,
                                        const std::vector<double>& stdev,
                                        const TrainState& s,
                                        const std::string& outFile,
                                        const std::string& moduleName)
{
    TFile fout(outFile.c_str(), "RECREATE");
    if (fout.IsZombie())
        throw cet::exception(moduleName) << "Cannot open partial-reverse diagnostic output file " << outFile;

    constexpr int kMaxDiagDims = 16;
    const int dim = static_cast<int>(data.front().x.size());
    if (dim > kMaxDiagDims)
        throw cet::exception(moduleName)
            << "Partial-reverse diagnostic supports at most " << kMaxDiagDims << " dimensions (got " << dim << ")";

    const size_t nUse = (s.partialReverseSamples > 0)
        ? std::min(data.size(), static_cast<size_t>(s.partialReverseSamples))
        : data.size();
    const size_t stride = std::max<size_t>(1, data.size() / nUse);
    const auto ranges = diagnosticTrueValueRanges(data, mean, stdev, dim, stride, nUse);
    constexpr int kNValBins = 200;
    constexpr int kN2DCap   = 500; // 2-D histograms are at most this many bins per axis
    const auto nbins = diagnosticBinCounts(dim, kNValBins, ranges, stdev, s.peakWindows);
    // Auto-range the residual histograms (created with empty limits) via the buffer mechanism;
    // their scale varies with t0. Restored before returning.
    const Int_t prevBufSize = TH1::GetDefaultBufferSize();
    TH1::SetDefaultBufferSize(50000);

    for (double t0 : s.partialReverseT0s) {
        std::ostringstream nm;
        nm << "prev_t" << t0;
        std::string treeName = nm.str();
        std::replace(treeName.begin(), treeName.end(), '.', 'p');
        std::ostringstream title;
        title << "Partial-reverse sampling diagnostic from t0=" << t0;
        TTree* tree = new TTree(treeName.c_str(), title.str().c_str());

        double t0Branch = t0;
        double trueD[kMaxDiagDims]    = {0.0};
        double sampledD[kMaxDiagDims] = {0.0};
        double noisedD[kMaxDiagDims]  = {0.0};
        tree->Branch("t0", &t0Branch);
        for (int i = 0; i < dim; ++i) {
            tree->Branch(("true_dim" + std::to_string(i)).c_str(), &trueD[i]);
            tree->Branch(("sampled_dim" + std::to_string(i)).c_str(), &sampledD[i]);
            tree->Branch(("noised_dim" + std::to_string(i)).c_str(), &noisedD[i]);
        }
        // Per-dim histograms: truth / noised start / sampled overlay (shared binning, finer on peak
        // dims), residual (sampled-truth, auto-ranged), and 2-D truth-vs-sampled (diagonal = perfect).
        std::vector<TH1D*> hTrue(dim), hSamp(dim), hNoised(dim), hResid(dim);
        std::vector<TH2D*> h2(dim);
        const std::string ts = std::to_string(t0);
        for (int i = 0; i < dim; ++i) {
            const std::string di = std::to_string(i);
            const int n2d = std::min(kN2DCap, nbins[i]);
            hTrue[i] = new TH1D((treeName + "_h_true_dim" + di).c_str(),
                ("truth dim" + di + " at t0=" + ts + ";value;count").c_str(),
                nbins[i], ranges[i].first, ranges[i].second);
            hSamp[i] = new TH1D((treeName + "_h_sampled_dim" + di).c_str(),
                ("partial-reverse sampled dim" + di + " at t0=" + ts + ";value;count").c_str(),
                nbins[i], ranges[i].first, ranges[i].second);
            hNoised[i] = new TH1D((treeName + "_h_noised_dim" + di).c_str(),
                ("noised start dim" + di + " at t0=" + ts + ";value;count").c_str(),
                nbins[i], ranges[i].first, ranges[i].second);
            hResid[i] = new TH1D((treeName + "_h_residual_dim" + di).c_str(),
                ("residual (sampled-truth) dim" + di + " at t0=" + ts + ";sampled - truth;count").c_str(),
                nbins[i], 0.0, 0.0); // auto-ranged via buffer
            h2[i] = new TH2D((treeName + "_h2_true_vs_sampled_dim" + di).c_str(),
                ("truth vs sampled dim" + di + " at t0=" + ts + ";truth;sampled").c_str(),
                n2d, ranges[i].first, ranges[i].second, n2d, ranges[i].first, ranges[i].second);
        }

        std::vector<double> noisedZ;
        size_t written = 0;
        for (size_t k = 0; k < data.size() && written < nUse; k += stride, ++written) {
            const auto& sample = data[k];
            auto res = model.partialReverseSample(sample.x, sample.cond, t0,
                s.denoiseDiagnosticUseEMA, s.partialReverseUseHeun, s.partialReverseUseSDE,
                s.partialReverseDiffusionSteps, s.partialReverseSdeToOdeSigmaThreshold, &noisedZ);
            for (int i = 0; i < dim; ++i) {
                trueD[i]    = sample.x[i] * stdev[i] + mean[i];
                sampledD[i] = res.value[i];
                noisedD[i]  = noisedZ[i] * stdev[i] + mean[i];
                hTrue[i]->Fill(trueD[i]);
                hSamp[i]->Fill(sampledD[i]);
                hNoised[i]->Fill(noisedD[i]);
                hResid[i]->Fill(sampledD[i] - trueD[i]);
                h2[i]->Fill(trueD[i], sampledD[i]);
            }
            tree->Fill();
        }
        tree->Write();
        for (int i = 0; i < dim; ++i) {
            hTrue[i]->Write(); hSamp[i]->Write(); hNoised[i]->Write();
            hResid[i]->BufferEmpty(); // flush buffer so the auto-range is computed before writing
            hResid[i]->Write();
            h2[i]->Write();
        }
        mf::LogInfo(moduleName)
            << "Partial-reverse diagnostic tree " << treeName << " filled (" << written << " samples)";
    }
    TH1::SetDefaultBufferSize(prevBufSize);
    fout.Close();
    mf::LogInfo(moduleName)
        << "Partial-reverse diagnostic written to " << outFile
        << " (" << nUse << " samples per t0, " << s.partialReverseT0s.size() << " t0 values)";
}

// ---------------------------------------------------------------------------
// runConditionalLossDiagnostic — profiles per-dimension loss vs sigma, split by whether each sample
//   lies inside a peak window, under TWO lenses: the eps-style loss (eps_hat-eps)^2 (all modes) and
//   the model's NATIVE training-target loss (EPS->eps, V->v, SCORE->score). Writes one TTree
//   "cond_loss" with per-sample rows: t, sigma, log10Sigma, windowIndex (the first peak window the
//   sample falls in, or -1), lossTotal + loss_dim<i> (eps-style), and nativeLossTotal +
//   nativeLoss_dim<i> (native target). Also writes, for each lens, TProfiles prof_LP_*/prof_LQ_* and
//   prof_LP_native_*/prof_LQ_native_* (loss vs log10 sigma) and ready-to-view overlay TCanvases:
//   c_LPLQ_* / c_LPLQ_native_* (L_P red vs L_Q blue) and c_epsVsNative_* (eps red vs native green,
//   out-of-window L_Q). Read L_P(sigma) vs L_Q(sigma): an underfit (L_P >> L_Q) only at low sigma
//   indicates a sampling/under-weighting problem (tune gMax/alpha); at all sigma it points to model
//   capacity / Fourier resolution. Under V the eps-style loss saturates ~1 at low sigma (eps is
//   unobservable there) so judge low-sigma fit by the native v-loss; the SCORE native loss scales
//   ~1/sigma^2 at low sigma and is not magnitude-comparable across sigma.
// ---------------------------------------------------------------------------
inline void runConditionalLossDiagnostic(ScoreBasedDiffusionModel& model,
                                         const std::vector<DiffusionTrainingSample>& data, // normalized
                                         const TrainState& s,
                                         const std::string& outFile,
                                         const std::string& moduleName)
{
    TFile fout(outFile.c_str(), "RECREATE");
    if (fout.IsZombie())
        throw cet::exception(moduleName) << "Cannot open conditional-loss diagnostic output file " << outFile;

    constexpr int kMaxDiagDims = 16;
    const int dim = static_cast<int>(data.front().x.size());
    if (dim > kMaxDiagDims)
        throw cet::exception(moduleName)
            << "Conditional-loss diagnostic supports at most " << kMaxDiagDims << " dimensions (got " << dim << ")";
    if (s.peakWindows.empty())
        mf::LogWarning(moduleName)
            << "Conditional-loss diagnostic: no peak windows configured; every row gets windowIndex=-1 "
            << "(L(sigma) over all data, no in/out-of-window split).";

    TTree* tree = new TTree("cond_loss", "Conditional loss vs sigma (eps-style + native target), split by feature window");
    double t = 0.0, sigma = 0.0, log10Sigma = 0.0, lossTotal = 0.0, nativeLossTotal = 0.0;
    int windowIndex = -1;
    double lossD[kMaxDiagDims] = {0.0};
    double nativeLossD[kMaxDiagDims] = {0.0};
    tree->Branch("t", &t);
    tree->Branch("sigma", &sigma);
    tree->Branch("log10Sigma", &log10Sigma);
    tree->Branch("windowIndex", &windowIndex);
    tree->Branch("lossTotal", &lossTotal);             // eps-style mean over dims
    tree->Branch("nativeLossTotal", &nativeLossTotal); // native-target mean over dims
    for (int i = 0; i < dim; ++i)
        tree->Branch(("loss_dim" + std::to_string(i)).c_str(), &lossD[i]);
    for (int i = 0; i < dim; ++i)
        tree->Branch(("nativeLoss_dim" + std::to_string(i)).c_str(), &nativeLossD[i]);

    // TProfiles of loss vs log10(sigma), split into in-window (P) and out-of-window (Q). The
    // log10-sigma axis spans the schedule's sigma range (t clamped away from the endpoints, as in
    // evalEpsLossSample).
    constexpr int kNSigBins = 60;
    const double logSigLo = std::log10(model.diffusionSigma(1e-3));
    const double logSigHi = std::log10(model.diffusionSigma(1.0 - 1e-3));
    auto makeProf = [&](const std::string& nm, const std::string& ti) {
        return new TProfile(nm.c_str(), ti.c_str(), kNSigBins, logSigLo, logSigHi);
    };
    TProfile* profLPtotal = makeProf("prof_LP_total", "L_P (in-window) mean eps-loss vs log_{10}#sigma;log_{10} #sigma;mean loss");
    TProfile* profLQtotal = makeProf("prof_LQ_total", "L_Q (out-of-window) mean eps-loss vs log_{10}#sigma;log_{10} #sigma;mean loss");
    std::vector<TProfile*> profLPdim(dim), profLQdim(dim);
    for (int i = 0; i < dim; ++i) {
        profLPdim[i] = makeProf("prof_LP_dim" + std::to_string(i),
            "L_P dim" + std::to_string(i) + " eps-loss vs log_{10}#sigma;log_{10} #sigma;mean loss");
        profLQdim[i] = makeProf("prof_LQ_dim" + std::to_string(i),
            "L_Q dim" + std::to_string(i) + " eps-loss vs log_{10}#sigma;log_{10} #sigma;mean loss");
    }
    // Native-target loss profiles (mode-dependent: EPS->eps, V->v, SCORE->score), parallel to the
    // eps-style ones above so the two lenses can be compared in the overlays below.
    TProfile* profLPnatTotal = makeProf("prof_LP_native_total", "L_P (in-window) mean native loss vs log_{10}#sigma;log_{10} #sigma;mean loss");
    TProfile* profLQnatTotal = makeProf("prof_LQ_native_total", "L_Q (out-of-window) mean native loss vs log_{10}#sigma;log_{10} #sigma;mean loss");
    std::vector<TProfile*> profLPnatDim(dim), profLQnatDim(dim);
    for (int i = 0; i < dim; ++i) {
        profLPnatDim[i] = makeProf("prof_LP_native_dim" + std::to_string(i),
            "L_P dim" + std::to_string(i) + " native loss vs log_{10}#sigma;log_{10} #sigma;mean loss");
        profLQnatDim[i] = makeProf("prof_LQ_native_dim" + std::to_string(i),
            "L_Q dim" + std::to_string(i) + " native loss vs log_{10}#sigma;log_{10} #sigma;mean loss");
    }

    // Window edges are in transformed (pre-z-score) units; convert to z-score once (the data is
    // z-scored) using the model's stored per-dimension mean/stdev.
    std::vector<double> zLow(s.peakWindows.size()), zHigh(s.peakWindows.size());
    for (size_t w = 0; w < s.peakWindows.size(); ++w) {
        zLow[w]  = model.normalizeCoord(s.peakWindows[w].dim, s.peakWindows[w].low);
        zHigh[w] = model.normalizeCoord(s.peakWindows[w].dim, s.peakWindows[w].high);
    }

    const size_t nUse = (s.condLossDiagnosticSamples > 0)
        ? std::min(data.size(), static_cast<size_t>(s.condLossDiagnosticSamples))
        : data.size();
    const size_t stride = std::max<size_t>(1, data.size() / nUse);
    size_t written = 0;
    for (size_t k = 0; k < data.size() && written < nUse; k += stride, ++written) {
        const auto& sample = data[k];
        windowIndex = -1;
        for (size_t w = 0; w < s.peakWindows.size(); ++w) {
            double v = sample.x[s.peakWindows[w].dim];
            if (v >= zLow[w] && v < zHigh[w]) { windowIndex = static_cast<int>(w); break; }
        }
        auto res = model.evalEpsLossSample(sample.x, sample.cond, s.denoiseDiagnosticUseEMA);
        t = res.t; sigma = res.sigma; log10Sigma = std::log10(res.sigma);
        lossTotal = 0.0; nativeLossTotal = 0.0;
        for (int i = 0; i < dim; ++i) {
            lossD[i]       = res.perDimLoss[i];       lossTotal       += res.perDimLoss[i];
            nativeLossD[i] = res.perDimNativeLoss[i]; nativeLossTotal += res.perDimNativeLoss[i];
        }
        lossTotal /= dim; nativeLossTotal /= dim;
        tree->Fill();

        const bool inWindow = (windowIndex >= 0);
        (inWindow ? profLPtotal : profLQtotal)->Fill(log10Sigma, lossTotal);
        (inWindow ? profLPnatTotal : profLQnatTotal)->Fill(log10Sigma, nativeLossTotal);
        for (int i = 0; i < dim; ++i) {
            (inWindow ? profLPdim[i]    : profLQdim[i])   ->Fill(log10Sigma, lossD[i]);
            (inWindow ? profLPnatDim[i] : profLQnatDim[i])->Fill(log10Sigma, nativeLossD[i]);
        }
    }
    tree->Write();
    profLPtotal->Write();
    profLQtotal->Write();
    for (int i = 0; i < dim; ++i) { profLPdim[i]->Write(); profLQdim[i]->Write(); }
    profLPnatTotal->Write();
    profLQnatTotal->Write();
    for (int i = 0; i < dim; ++i) { profLPnatDim[i]->Write(); profLQnatDim[i]->Write(); }

    // L_P vs L_Q overlay canvas per pair (total and each dim), so the comparison is visible on
    // opening the file. Force batch mode while building canvases so this is safe on headless nodes.
    const bool prevBatch = gROOT->IsBatch();
    gROOT->SetBatch(kTRUE);
    auto writeOverlay = [&](const std::string& nm, TProfile* lp, TProfile* lq) {
        lp->SetStats(0); lq->SetStats(0); // no statbox on the overlay
        lp->SetLineColor(kRed);  lp->SetMarkerColor(kRed);  lp->SetMarkerStyle(20);
        lq->SetLineColor(kBlue); lq->SetMarkerColor(kBlue); lq->SetMarkerStyle(21);
        TCanvas c(nm.c_str(), nm.c_str(), 800, 600);
        const double ymax = std::max(lp->GetMaximum(), lq->GetMaximum());
        lp->SetMinimum(0.0);
        lp->SetMaximum(ymax > 0.0 ? 1.1 * ymax : 1.0);
        lp->Draw();
        lq->Draw("SAME");
        TLegend leg(0.66, 0.76, 0.88, 0.88);
        leg.AddEntry(lp, "L_P (in-window)",     "lep");
        leg.AddEntry(lq, "L_Q (out-of-window)", "lep");
        leg.Draw();
        c.Write();
    };
    writeOverlay("c_LPLQ_total", profLPtotal, profLQtotal);
    for (int i = 0; i < dim; ++i)
        writeOverlay("c_LPLQ_dim" + std::to_string(i), profLPdim[i], profLQdim[i]);
    // Native-target L_P vs L_Q overlays, paralleling the eps-style ones above.
    writeOverlay("c_LPLQ_native_total", profLPnatTotal, profLQnatTotal);
    for (int i = 0; i < dim; ++i)
        writeOverlay("c_LPLQ_native_dim" + std::to_string(i), profLPnatDim[i], profLQnatDim[i]);

    // Eps-style vs native-target overlay (out-of-window L_Q lens), so the two loss definitions sit on
    // one canvas. Under V the eps curve saturates ~1 at low sigma while native (v) stays informative;
    // under EPS the two coincide. CLONE both profiles before recoloring: the originals are already
    // referenced by the LP/LQ canvases written above, so mutating their draw attributes here would
    // retroactively corrupt those canvases on read-back. Clones are owned by (and freed with) the
    // canvas via kCanDelete.
    auto writeEpsVsNative = [&](const std::string& nm, TProfile* epsProf, TProfile* natProf) {
        TCanvas c(nm.c_str(), nm.c_str(), 800, 600);
        auto* epsC = static_cast<TProfile*>(epsProf->Clone((nm + "_eps").c_str()));
        auto* natC = static_cast<TProfile*>(natProf->Clone((nm + "_native").c_str()));
        // Detach from the open TFile (Clone auto-registers) so the canvas is the sole owner; then
        // hand ownership to the canvas via kCanDelete. Avoids a double-free and stray file objects.
        epsC->SetDirectory(nullptr); natC->SetDirectory(nullptr);
        epsC->SetBit(kCanDelete); natC->SetBit(kCanDelete);
        epsC->SetStats(0); natC->SetStats(0);
        epsC->SetLineColor(kRed);   epsC->SetMarkerColor(kRed);   epsC->SetMarkerStyle(20);
        natC->SetLineColor(kGreen + 2); natC->SetMarkerColor(kGreen + 2); natC->SetMarkerStyle(22);
        const double ymax = std::max(epsC->GetMaximum(), natC->GetMaximum());
        epsC->SetMinimum(0.0);
        epsC->SetMaximum(ymax > 0.0 ? 1.1 * ymax : 1.0);
        epsC->Draw();
        natC->Draw("SAME");
        TLegend leg(0.60, 0.76, 0.88, 0.88);
        leg.AddEntry(epsC, "eps-style loss",   "lep");
        leg.AddEntry(natC, "native-target loss", "lep");
        leg.Draw();
        c.Write();
    };
    writeEpsVsNative("c_epsVsNative_total", profLQtotal, profLQnatTotal);
    for (int i = 0; i < dim; ++i)
        writeEpsVsNative("c_epsVsNative_dim" + std::to_string(i), profLQdim[i], profLQnatDim[i]);
    gROOT->SetBatch(prevBatch);

    fout.Close();
    mf::LogInfo(moduleName)
        << "Conditional-loss diagnostic written to " << outFile << " (" << written << " samples, "
        << s.peakWindows.size() << " window(s)).";
}

// ---------------------------------------------------------------------------
// runFeatureBlockDiagnostic — first-layer input-feature-block weight/grad L2 magnitudes plus a
//   fresh per-output-dimension loss. Writes TTree "feature_blocks" (one row per block: blockIndex,
//   kind, coord, nCols, weightL2, gradL2, name) and TTree "output_dim_loss" (dim, loss), and logs a
//   readable table. Near-zero gradL2 / init-scale weightL2 on a block (e.g. fourier_state[pz]) means
//   that feature is not being used or driven — distinguishing a representation problem from a
//   sampling one.
// ---------------------------------------------------------------------------
inline void runFeatureBlockDiagnostic(ScoreBasedDiffusionModel& model,
                                      const std::vector<DiffusionTrainingSample>& data, // normalized
                                      const TrainState& s,
                                      const std::string& outFile,
                                      const std::string& moduleName)
{
    std::vector<double> perDimLoss;
    auto blocks = model.firstLayerBlockMagnitudes(data, s.featureBlockDiagnosticSamples, perDimLoss);

    TFile fout(outFile.c_str(), "RECREATE");
    if (fout.IsZombie())
        throw cet::exception(moduleName) << "Cannot open feature-block diagnostic output file " << outFile;

    TTree* bt = new TTree("feature_blocks", "First-layer input-feature-block weight/grad L2");
    int blockIndex = 0, kind = 0, coord = 0, nCols = 0;
    double weightL2 = 0.0, gradL2 = 0.0;
    std::string name;
    bt->Branch("blockIndex", &blockIndex);
    bt->Branch("kind", &kind);
    bt->Branch("coord", &coord);
    bt->Branch("nCols", &nCols);
    bt->Branch("weightL2", &weightL2);
    bt->Branch("gradL2", &gradL2);
    bt->Branch("name", &name);

    // Per-block bar charts (bin label = block name) so the file is self-contained.
    const int nb = static_cast<int>(blocks.size());
    TH1D* hWeight = new TH1D("h_weightL2", "First-layer weight L2 per input-feature block;;weight L2", nb, 0, nb);
    TH1D* hGrad   = new TH1D("h_gradL2",   "First-layer mean-gradient L2 per input-feature block;;grad L2", nb, 0, nb);
    hWeight->SetStats(0); hGrad->SetStats(0); // no statbox on the bar charts
    std::ostringstream tbl;
    tbl << "First-layer feature-block magnitudes "
        << "(kind: 0 raw-state 1 fourier-state 2 raw-cond 3 fourier-cond 4 raw-time 5 fourier-time):";
    for (size_t b = 0; b < blocks.size(); ++b) {
        blockIndex = static_cast<int>(b);
        kind = blocks[b].kind; coord = blocks[b].coord; nCols = blocks[b].nCols;
        weightL2 = blocks[b].weightL2; gradL2 = blocks[b].gradL2; name = blocks[b].name;
        bt->Fill();
        hWeight->SetBinContent(static_cast<int>(b) + 1, blocks[b].weightL2);
        hWeight->GetXaxis()->SetBinLabel(static_cast<int>(b) + 1, blocks[b].name.c_str());
        hGrad->SetBinContent(static_cast<int>(b) + 1, blocks[b].gradL2);
        hGrad->GetXaxis()->SetBinLabel(static_cast<int>(b) + 1, blocks[b].name.c_str());
        tbl << "\n  " << name << "  nCols=" << nCols << "  weightL2=" << weightL2 << "  gradL2=" << gradL2;
    }
    bt->Write();
    hWeight->Write();
    hGrad->Write();
    mf::LogInfo(moduleName) << tbl.str();

    TTree* dt = new TTree("output_dim_loss", "Fresh mean per-output-dimension squared residual");
    int dimIndex = 0; double dimLoss = 0.0;
    dt->Branch("dim", &dimIndex);
    dt->Branch("loss", &dimLoss);
    const int nd = static_cast<int>(perDimLoss.size());
    TH1D* hDimLoss = new TH1D("h_output_dim_loss", "Per-output-dimension mean squared residual;dimension;mean loss", nd, 0, nd);
    hDimLoss->SetStats(0); // no statbox
    std::ostringstream dloss; dloss << "Per-output-dimension mean squared residual:";
    for (size_t i = 0; i < perDimLoss.size(); ++i) {
        dimIndex = static_cast<int>(i); dimLoss = perDimLoss[i];
        dt->Fill();
        hDimLoss->SetBinContent(static_cast<int>(i) + 1, perDimLoss[i]);
        dloss << "\n  dim" << i << " = " << perDimLoss[i];
    }
    dt->Write();
    hDimLoss->Write();
    fout.Close();
    mf::LogInfo(moduleName) << dloss.str();
    mf::LogInfo(moduleName) << "Feature-block diagnostic written to " << outFile;
}

// ---------------------------------------------------------------------------
// runTraining — finalise normalisation, then run the full per-epoch training
//               loop (with curriculum phase transitions and epoch checkpointing)
//               for each active model path. If denoiseDiagnosticTs is set, the
//               one-step denoising diagnostic replaces training entirely.
// ---------------------------------------------------------------------------
inline void runTraining(TrainState& s, const std::string& moduleName) {
    // Finalise Welford stdevs
    auto stdev = [&](double M2) { return std::sqrt(M2 / static_cast<double>(s.nNorm)); };
    s.t_stdev    = stdev(s.t_M2);
    s.x_stdev    = stdev(s.x_M2);
    s.y_stdev    = stdev(s.y_M2);
    s.mom0_stdev = stdev(s.mom0_M2);
    s.mom1_stdev = stdev(s.mom1_M2);
    s.mom2_stdev = stdev(s.mom2_M2);

    // Guard against empty datasets
    if (s.useTwoStageTraining) {
        if (s.trainStage1 && s.stage1TrainingData.empty()) {
            mf::LogWarning(moduleName) << "No training data collected for stage-1.";
            return;
        }
        if (s.trainStage2 && s.stage2TrainingData.empty()) {
            mf::LogWarning(moduleName) << "No training data collected for stage-2.";
            return;
        }
    } else if (s.allAtOnceTrainingData.empty()) {
        mf::LogWarning(moduleName) << "No training data collected.";
        return;
    }

    // Strip trailing ".csv" or ".bin" for checkpoint naming
    auto stripExt = [](const std::string& f) -> std::string {
        if (f.size() > 4 && (f.substr(f.size() - 4) == ".csv" || f.substr(f.size() - 4) == ".bin"))
            return f.substr(0, f.size() - 4);
        return f;
    };

    // Diagnostic mode: when either diagnostic is requested it REPLACES training for
    // every active model path. Shared by all three paths below.
    const bool diagnosticMode = !s.denoiseDiagnosticTs.empty() || !s.partialReverseT0s.empty()
                                || s.condLossDiagnostic || s.featureBlockDiagnostic;
    auto runDiagnostics = [&](ScoreBasedDiffusionModel& model,
                              const std::vector<DiffusionTrainingSample>& data,
                              const std::vector<double>& mean,
                              const std::vector<double>& stdevs,
                              const std::string& modelFile) {
        if (!s.denoiseDiagnosticTs.empty())
            runDenoiseDiagnostic(model, data, mean, stdevs, s,
                stripExt(modelFile) + "_DenoisingDiag.root", moduleName);
        if (!s.partialReverseT0s.empty())
            runPartialReverseDiagnostic(model, data, mean, stdevs, s,
                stripExt(modelFile) + "_PartialReverseDiag.root", moduleName);
        if (s.condLossDiagnostic)
            runConditionalLossDiagnostic(model, data, s,
                stripExt(modelFile) + "_CondLossDiag.root", moduleName);
        if (s.featureBlockDiagnostic)
            runFeatureBlockDiagnostic(model, data, s,
                stripExt(modelFile) + "_FeatureBlockDiag.root", moduleName);
    };

    // Per-epoch training loop shared by all three model paths.
    // resetPhase0: if true, resets currentBias/tLow to phase-0 values before the loop
    //              (needed for stage-2, which follows stage-1 in the same endJob call).
    auto trainLoop = [&](ScoreBasedDiffusionModel& model,
                         std::vector<DiffusionTrainingSample>& data,
                         const std::string& outFile,
                         bool resetPhase0,
                         int basisTag) { // opaque tag persisted with every saved checkpoint
        // Set every peak window's alpha for curriculum phase k: the per-phase override
        // (curriculumPeakAlpha[k]) when >= 0, else the window's snapshotted base alpha.
        // train() re-reads s.peakWindows each epoch, so mutating it here takes effect for
        // the phase. Returns the effective alpha for logging (or -1 if no windows). No-op
        // under a non-curriculum run (nPhase <= 1) where windows keep their base alpha.
        auto applyPeakAlpha = [&](int k) -> double {
            if (s.nPhase <= 1 || s.peakWindows.empty()) return -1.0;
            const double aOverride = s.curriculumPeakAlpha[k];
            double eff = -1.0;
            for (size_t w = 0; w < s.peakWindows.size(); ++w) {
                s.peakWindows[w].alpha = (aOverride >= 0.0) ? aOverride : s.basePeakAlphas[w];
                eff = s.peakWindows[w].alpha;
            }
            return eff;
        };
        // Set the live per-phase window list fed to train(): the configured windows when phase k
        // enables peak sampling, else empty (uniform draws, wIS=1). A non-curriculum run keeps
        // the windows on. peakAlpha for a sampling-off phase becomes moot (no in-window draws).
        auto applyPeakSampling = [&](int k) {
            const bool on = (s.nPhase <= 1) ? true : s.curriculumUsePeakSampling[k];
            s.currentPeakWindows = on ? s.peakWindows : std::vector<PeakWindow>{};
        };
        if (resetPhase0) {
            s.currentBiasLowSigma    = s.nPhase > 1 ? s.curriculumBiasLowSigma[0]    : s.currentBiasLowSigma;
            s.currentTLowBound       = s.nPhase > 1 ? s.curriculumTLowBound[0]       : s.currentTLowBound;
            s.currentTFocusLow       = s.nPhase > 1 ? s.curriculumTFocusLow[0]       : s.currentTFocusLow;
            s.currentTFocusHigh      = s.nPhase > 1 ? s.curriculumTFocusHigh[0]      : s.currentTFocusHigh;
            s.currentTFocusFraction  = s.nPhase > 1 ? s.curriculumTFocusFraction[0]  : s.currentTFocusFraction;
            s.currentSubsetSizePerEpoch = s.nPhase > 1 ? s.curriculumSubsetSizePerEpoch[0] : s.currentSubsetSizePerEpoch;
            s.currentUsePeakWindowLoss  = s.nPhase > 1 ? s.curriculumUsePeakWindowLoss[0]  : s.currentUsePeakWindowLoss;
        }
        // Apply phase-0 peak alpha unconditionally (override or base): enterPhase is only
        // called for k>=1, and resetPhase0 may be false on the primary all-at-once/stage-1
        // paths, so phase 0's override would otherwise never take effect.
        applyPeakAlpha(0);
        applyPeakSampling(0); // phase-0 peak-sampling on/off (also seeds currentPeakWindows)
        if (s.promoteEMAOnStart || (s.nPhase > 1 && s.curriculumPromoteEMA[0]))
            model.promoteEMAToNetwork();
        mf::LogInfo(moduleName)
            << "Initial training settings: biasLowSigma=" << s.currentBiasLowSigma
            << ", tLowBound=" << s.currentTLowBound
            << ", tFocusWindow=[" << s.currentTFocusLow << ", " << s.currentTFocusHigh
            << "] (fraction=" << s.currentTFocusFraction << ")"
            << ", trainingSubsetSizePerEpoch=" << s.currentSubsetSizePerEpoch;
        const std::string base = stripExt(outFile);

        // Apply curriculum phase `k`'s hyperparameters (k is 0-based; only meaningful
        // for k>=1 — phase 0 keeps the constructed-model defaults, as in legacy mode).
        // Logs "Switching to phase k+1" to match the existing 1-based phase numbering.
        auto enterPhase = [&](int k) {
            double lwp = s.curriculumLossWeightPower[k];
            double gc  = s.curriculumGradientClip   [k];
            double lr  = s.curriculumLearningRate   [k];
            int    bs  = s.curriculumBatchSize      [k];
            bool   udwc= s.curriculumUseDimWeightController[k];
            s.currentBiasLowSigma    = s.curriculumBiasLowSigma   [k];
            s.currentTLowBound       = s.curriculumTLowBound      [k];
            s.currentTFocusLow       = s.curriculumTFocusLow      [k];
            s.currentTFocusHigh      = s.curriculumTFocusHigh     [k];
            s.currentTFocusFraction  = s.curriculumTFocusFraction [k];
            s.currentSubsetSizePerEpoch = s.curriculumSubsetSizePerEpoch[k];
            s.currentUsePeakWindowLoss  = s.curriculumUsePeakWindowLoss[k];
            const double effAlpha = applyPeakAlpha(k); // mutate s.peakWindows alpha for this phase
            applyPeakSampling(k); // set currentPeakWindows (windows or empty) for this phase
            mf::LogInfo(moduleName)
                << "Switching to phase " << (k + 1)
                << ": lossWeightPower=" << lwp
                << ", gradientClipThreshold=" << gc
                << ", learningRate=" << lr
                << ", batchSize=" << bs
                << ", useDimWeightController=" << udwc
                << ", biasLowSigma=" << s.currentBiasLowSigma
                << ", tLowBound=" << s.currentTLowBound
                << ", tFocusWindow=[" << s.currentTFocusLow << ", " << s.currentTFocusHigh
                << "] (fraction=" << s.currentTFocusFraction << ")"
                << ", trainingSubsetSizePerEpoch=" << s.currentSubsetSizePerEpoch
                << ", usePeakWindowLoss=" << s.currentUsePeakWindowLoss
                << ", peakAlpha=" << (s.curriculumPeakAlpha[k] >= 0.0 ? effAlpha : -1.0)
                << (s.curriculumPeakAlpha[k] >= 0.0 ? " (override)" : " (base)")
                << ", peakSampling=" << (s.currentPeakWindows.empty() ? "off" : "on");
            model.updateLossWeightPower(lwp);
            model.updateGradientClipThreshold(gc);
            model.updateLearningRate(lr);
            model.updateBatchSize(bs);
            // Apply before promoteEMA so the frozen-weights log (if the controller
            // turns off here) reflects the pre-promotion state.
            model.updateUseDimWeightController(udwc);
            if (s.curriculumPromoteEMA[k])
                model.promoteEMAToNetwork();
        };

        if (!s.autoPlanner) {
            // ----- Legacy fixed-epoch schedule (unchanged behavior) -----
            for (int e = 1; e <= s.trainingEpochs; ++e) {
                if (s.nPhase > 1) {
                    for (int ph = 0; ph < s.nPhase - 1; ++ph)
                        if (e == s.phaseBoundaries[ph] + 1) enterPhase(ph + 1);
                }
                mf::LogInfo(moduleName) << "Epoch " << e << "/" << s.trainingEpochs;
                model.train(data, 1, s.currentSubsetSizePerEpoch, s.currentBiasLowSigma, s.currentTLowBound,
                            s.currentTFocusLow, s.currentTFocusHigh, s.currentTFocusFraction, s.currentPeakWindows);
                if (std::find(s.saveEpochs.begin(), s.saveEpochs.end(), e) != s.saveEpochs.end()) {
                    model.saveModel(base + ".epoch" + std::to_string(e) + ".bin", basisTag);
                    if (s.saveAlsoCsv)
                        model.saveModelCsv(base + ".epoch" + std::to_string(e) + ".csv");
                }
            }
            model.saveModel(outFile, basisTag);
            if (s.saveAlsoCsv)
                model.saveModelCsv(base + ".csv");
            mf::LogInfo(moduleName) << "Training completed, saved to " << outFile;
            return;
        }

        // ----- Automatic curriculum planner: train each phase until converged -----
        std::string capDesc;
        if (s.nPhase > 1) {
            std::ostringstream oss;
            oss << "SBDMtrainingCurriculumEpochs=[";
            for (int i = 0; i < s.nPhase; ++i) oss << (i ? ", " : "") << s.curriculumEpochs[i];
            oss << "] (0 = skipped phase)";
            capDesc = oss.str();
        } else {
            capDesc = "SBDMplannerMaxEpochsPerPhase=" + std::to_string(s.plannerMaxEpochsPerPhase);
        }
        std::string minDeltaDesc;
        if (s.nPhase > 1) {
            std::ostringstream oss;
            oss << "SBDMtrainingCurriculumMinDelta=[";
            for (int i = 0; i < s.nPhase; ++i) oss << (i ? ", " : "") << s.curriculumMinDelta[i];
            oss << "]";
            minDeltaDesc = oss.str();
        } else {
            minDeltaDesc = std::to_string(s.plannerMinDelta);
        }
        mf::LogInfo(moduleName)
            << "Auto curriculum planner enabled: smoothWindow=" << s.plannerSmoothWindow
            << ", minDelta=" << minDeltaDesc << ", patience=" << s.plannerPatience
            << ", minEpochsPerPhase=" << s.plannerMinEpochsPerPhase
            << ", max-epoch cap per phase=" << capDesc;
        ConvergenceTracer tracer{s.plannerSmoothWindow, s.plannerMinDelta,
                                   s.plannerPatience, s.plannerMinEpochsPerPhase};
        int globalEpoch = 0;
        for (int ph = 0; ph < s.nPhase; ++ph) {
            const int phaseNum = ph + 1; // base-1 for logs and filenames
            // A phase configured with 0 epochs is a "setup-only" phase: apply its
            // hyperparameters (for ph>=1) and advance without training. This preserves
            // the resume idiom SBDMtrainingCurriculumEpochs=[0, N], where phase 0 exists
            // only to push the real phase-1 settings (lr, clip, promoteEMA, ...) onto a
            // model loaded from an old checkpoint before any new training happens.
            if (s.nPhase > 1 && s.curriculumEpochs[ph] == 0) {
                if (ph >= 1) enterPhase(ph);
                mf::LogInfo(moduleName)
                    << "Phase " << phaseNum << " configured with 0 epochs; applied its "
                    << "settings and skipped (no training).";
                continue;
            }
            if (ph >= 1) enterPhase(ph);
            tracer.reset();
            tracer.minDelta = (s.nPhase > 1) ? s.curriculumMinDelta[ph] : s.plannerMinDelta;
            model.clearNetworkSnapshot(); // best-of-phase snapshot is per-phase
            // Per-phase max-epoch cap: a curriculum phase's configured epoch count
            // (guaranteed > 0 here, since 0-epoch phases were skipped above); for a
            // single-phase / no-curriculum run, the global planner cap. Validated in
            // validateAndBuildCurriculum to be >= the convergence floor.
            const int maxEpochs = (s.nPhase > 1) ? s.curriculumEpochs[ph]
                                                 : s.plannerMaxEpochsPerPhase;
            const std::string phaseTag = base + ".phase" + std::to_string(phaseNum);

            // Diagnostic: loss at the weights as just installed by the prior phase's
            // restoreNetwork() and/or this phase's enterPhase() promote — BEFORE any
            // optimizer step. A sane value here with an exploding epoch-1 training loss
            // means the weights are fine and the first step is to blame; a huge value
            // means the installed weights themselves are bad.
            {
                double evalLoss = model.evaluateAverageLoss(
                    data, s.currentSubsetSizePerEpoch, s.currentBiasLowSigma, s.currentTLowBound,
                    s.currentTFocusLow, s.currentTFocusHigh, s.currentTFocusFraction);
                mf::LogInfo(moduleName)
                    << "Phase " << phaseNum << " pre-train eval loss (no optimizer step) = " << evalLoss;
            }

            bool converged = false;
            int  epochInPhase = 0;
            while (true) {
                ++epochInPhase;
                ++globalEpoch;
                mf::LogInfo(moduleName)
                    << "Phase " << phaseNum << " epoch " << epochInPhase
                    << " (global epoch " << globalEpoch << ")";
                model.train(data, 1, s.currentSubsetSizePerEpoch, s.currentBiasLowSigma, s.currentTLowBound,
                            s.currentTFocusLow, s.currentTFocusHigh, s.currentTFocusFraction, s.currentPeakWindows);
                const double aggLoss = model.getLastEpochLoss(); // read immediately after train()
                // Metric selection: in a phase flagged usePeakWindowLoss, plateau on the
                // peak-window (feature-region) loss when it is available; fall back to the
                // aggregate loss on any epoch with no in-window draws (metric NaN).
                const double peakLoss  = model.getLastEpochPeakLoss();
                const long   peakCount = model.getLastEpochPeakCount();
                const bool   usePeak   = s.currentUsePeakWindowLoss
                                       && std::isfinite(peakLoss) && peakCount > 0;
                const double loss      = usePeak ? peakLoss : aggLoss; // scalar fed to the tracer
                const bool conv = tracer.update(loss);

                // Per-epoch planner diagnostics: smoothed loss, its relative improvement
                // vs the running best (the "delta" compared against minDelta), and the
                // no-improvement streak toward convergence. Makes convergence behavior
                // readable from the log instead of reverse-engineered from raw Loss.
                {
                    std::ostringstream pss;
                    pss << "Phase " << phaseNum << " planner: tracking="
                        << (usePeak ? "peakWindowLoss" : "aggregateLoss")
                        << " aggLoss=" << aggLoss
                        << " peakLoss=" << peakLoss << " peakCount=" << peakCount
                        << " rawLoss=" << loss
                        << " smoothed=" << tracer.lastSmoothed;
                    if (std::isfinite(tracer.lastDelta))
                        pss << " delta=" << (tracer.lastDelta * 100.0) << "% (minDelta="
                            << (tracer.minDelta * 100.0) << "%)";
                    else
                        pss << " delta=n/a";
                    pss << " sinceImprove=" << tracer.sinceImprove << "/" << s.plannerPatience
                        << " bestSmoothed=" << tracer.bestSmoothed
                        << " @epoch " << tracer.bestSmoothedEpoch;
                    mf::LogInfo(moduleName) << pss.str();
                }

                // Best-in-phase checkpoints. Raw = the single lowest-loss epoch;
                // smoothed = lowest moving-average loss (robust to per-epoch noise).
                if (tracer.rawImproved) {
                    model.saveModel(phaseTag + ".bestRaw.bin", basisTag);
                    if (s.saveAlsoCsv)
                        model.saveModelCsv(phaseTag + ".bestRaw.csv");
                }
                if (tracer.smoothedImproved) {
                    // Capture the smoothed-best in memory for resume-from-best, and on disk.
                    model.snapshotNetwork();
                    model.saveModel(phaseTag + ".bestSmoothed.bin", basisTag);
                    if (s.saveAlsoCsv)
                        model.saveModelCsv(phaseTag + ".bestSmoothed.csv");
                }
                // Honor explicit checkpoint epochs (indexed by global epoch count).
                if (std::find(s.saveEpochs.begin(), s.saveEpochs.end(), globalEpoch) != s.saveEpochs.end()) {
                    model.saveModel(base + ".epoch" + std::to_string(globalEpoch) + ".bin", basisTag);
                    if (s.saveAlsoCsv)
                        model.saveModelCsv(base + ".epoch" + std::to_string(globalEpoch) + ".csv");
                }

                if (conv) { converged = true; break; }
                if (epochInPhase >= maxEpochs) break;
            }

            // Phase-final snapshot: the literal final-epoch state, saved BEFORE any
            // resume-from-best restore so ".final" always means "last epoch trained".
            model.saveModel(phaseTag + ".final.bin", basisTag);
            if (s.saveAlsoCsv)
                model.saveModelCsv(phaseTag + ".final.csv");

            // Resume-from-best: ALWAYS rewind to this phase's smoothed-best before advancing.
            // The best point is the consistent reference for the whole training state — base
            // network + Adam moments, the EMA copy, and the dim-weight controller. Doing this
            // unconditionally (not just when the next phase skips EMA promotion) ensures that:
            //  - a non-promoting next phase / the global-final model continues from the best
            //    base weights (not the patience-drifted final epoch), and
            //  - a promoting next phase promotes the best-time EMA (restored here) with the
            //    best-time controller state, instead of the drifted phase-end EMA/dimWeights.
            if (model.hasNetworkSnapshot()) {
                model.restoreNetwork();
                const bool promoteNext = (ph + 1 < s.nPhase) ? s.curriculumPromoteEMA[ph + 1] : false;
                mf::LogInfo(moduleName)
                    << "Phase " << phaseNum << ": restored smoothed-best state (epoch "
                    << tracer.bestSmoothedEpoch << ", smoothed loss " << tracer.bestSmoothed << ")"
                    << (promoteNext ? " — next phase will promote the restored (best-time) EMA."
                                    : (ph + 1 < s.nPhase ? " as the next phase's starting point."
                                                         : " as the final model."));
            }

            auto bestSummary = [&]() {
                std::stringstream ss;
                ss << "raw-best epoch=" << tracer.bestRawEpoch << " (loss=" << tracer.bestRawLoss
                   << "), smoothed-best epoch=" << tracer.bestSmoothedEpoch
                   << " (smoothed loss=" << tracer.bestSmoothed << ")";
                return ss.str();
            };
            if (converged) {
                mf::LogInfo(moduleName)
                    << "Phase " << phaseNum << " converged after " << epochInPhase
                    << " epochs; " << bestSummary() << ".";
            } else {
                mf::LogWarning(moduleName)
                    << "Phase " << phaseNum << " did NOT converge within its max-epoch cap="
                    << maxEpochs << " epochs (" << bestSummary()
                    << "). Adjust the training plan — e.g. raise the cap, change the learning "
                    << "rate, or revisit the phase hyperparameters. Stopping training.";
                model.saveModel(outFile, basisTag);
                if (s.saveAlsoCsv)
                    model.saveModelCsv(base + ".csv");
                return; // abort: surface the bad plan instead of advancing
            }
        }
        model.saveModel(outFile, basisTag);
        if (s.saveAlsoCsv)
            model.saveModelCsv(base + ".csv");
        mf::LogInfo(moduleName) << "Training completed (auto planner), saved to " << outFile;
    };

    // Per-stage normalization mean/stdev orderings MUST match the vectors built in
    // collectSample for the active basis. V1 keeps the legacy arrangement; V2 puts
    // pTotal (mom0) alone in stage1 and conditions the 5-D stage2 on it.
    const bool isV2 = s.isV2();
    if (s.useTwoStageTraining) {
        if (s.trainStage1) {
            if (s.trainingSize > 0 && (int)s.stage1TrainingData.size() > s.trainingSize)
                s.stage1TrainingData.resize(s.trainingSize);
            // V2 stage1 = {mom0}; V1 stage1 = {t,x,y}.
            std::vector<double> mean  = isV2 ? std::vector<double>{s.mom0_mean}
                                             : std::vector<double>{s.t_mean, s.x_mean, s.y_mean};
            std::vector<double> sdev  = isV2 ? std::vector<double>{s.mom0_stdev}
                                             : std::vector<double>{s.t_stdev, s.x_stdev, s.y_stdev};
            s.stage1Model->normalizeData(mean, sdev, s.stage1TrainingData);
            const int tag1 = packBasisTag(isV2 ? ModelLayout::TwoStageStage1Ptot1D
                                               : ModelLayout::AllAtOnce6D, s.momentumBasis);
            if (diagnosticMode) {
                runDiagnostics(*s.stage1Model, s.stage1TrainingData, mean, sdev, s.stage1ModelFile);
            } else {
                mf::LogInfo(moduleName)
                    << "Training stage-1 diffusion model with " << s.stage1TrainingData.size()
                    << " samples and " << s.trainingEpochs << " epochs...";
                trainLoop(*s.stage1Model, s.stage1TrainingData, s.stage1ModelFile, false, tag1);
            }
        }
        if (s.trainStage2) {
            if (s.trainingSize > 0 && (int)s.stage2TrainingData.size() > s.trainingSize)
                s.stage2TrainingData.resize(s.trainingSize);
            // normalizeData order = [data dims..., cond dims...].
            // V2 stage2 data {t,x,y,mom1,mom2} cond {mom0};
            // V1 stage2 data {mom0,mom1,mom2} cond {t,x,y}.
            std::vector<double> mean = isV2
                ? std::vector<double>{s.t_mean, s.x_mean, s.y_mean, s.mom1_mean, s.mom2_mean, s.mom0_mean}
                : std::vector<double>{s.mom0_mean, s.mom1_mean, s.mom2_mean, s.t_mean, s.x_mean, s.y_mean};
            std::vector<double> sdev = isV2
                ? std::vector<double>{s.t_stdev, s.x_stdev, s.y_stdev, s.mom1_stdev, s.mom2_stdev, s.mom0_stdev}
                : std::vector<double>{s.mom0_stdev, s.mom1_stdev, s.mom2_stdev, s.t_stdev, s.x_stdev, s.y_stdev};
            s.stage2Model->normalizeData(mean, sdev, s.stage2TrainingData);
            const int tag2 = packBasisTag(isV2 ? ModelLayout::TwoStageStage2_5D
                                               : ModelLayout::AllAtOnce6D, s.momentumBasis);
            if (diagnosticMode) {
                runDiagnostics(*s.stage2Model, s.stage2TrainingData, mean, sdev, s.stage2ModelFile);
            } else {
                mf::LogInfo(moduleName)
                    << "Training stage-2 diffusion model with " << s.stage2TrainingData.size()
                    << " samples and " << s.trainingEpochs << " epochs...";
                trainLoop(*s.stage2Model, s.stage2TrainingData, s.stage2ModelFile, true, tag2);
            }
        }
    } else {
        if (s.allAtOnceModelFile.empty())
            throw cet::exception(moduleName) << "All-at-once training requires SBDMallAtOnceModelFile.";
        if (s.trainingSize > 0 && (int)s.allAtOnceTrainingData.size() > s.trainingSize)
            s.allAtOnceTrainingData.resize(s.trainingSize);
        // All-at-once vector (t,x,y,mom0,mom1,mom2) in either basis.
        std::vector<double> mean = {s.t_mean, s.x_mean, s.y_mean, s.mom0_mean, s.mom1_mean, s.mom2_mean};
        std::vector<double> sdev = {s.t_stdev, s.x_stdev, s.y_stdev, s.mom0_stdev, s.mom1_stdev, s.mom2_stdev};
        s.allAtOnceModel->normalizeData(mean, sdev, s.allAtOnceTrainingData);
        const int tagA = packBasisTag(ModelLayout::AllAtOnce6D, s.momentumBasis);
        if (diagnosticMode) {
            runDiagnostics(*s.allAtOnceModel, s.allAtOnceTrainingData, mean, sdev, s.allAtOnceModelFile);
        } else {
            mf::LogInfo(moduleName)
                << "Training all-at-once diffusion model with " << s.allAtOnceTrainingData.size()
                << " samples and " << s.trainingEpochs << " epochs...";
            trainLoop(*s.allAtOnceModel, s.allAtOnceTrainingData, s.allAtOnceModelFile, false, tagA);
        }
    }
}

} // namespace VDResampler
} // namespace mu2e
