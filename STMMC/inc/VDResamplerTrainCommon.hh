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
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TFile.h"
#include "TTree.h"

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

    // Welford normalization accumulators (transformed: t, x, y, pr, pphi, pz)
    double t_mean=0,    t_M2=0,    t_stdev=0;
    double x_mean=0,    x_M2=0,    x_stdev=0;
    double y_mean=0,    y_M2=0,    y_stdev=0;
    double pr_mean=0,   pr_M2=0,   pr_stdev=0;
    double pphi_mean=0, pphi_M2=0, pphi_stdev=0;
    double pz_mean=0,   pz_M2=0,   pz_stdev=0;
    int nNorm = 0;
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
    bool   epsPrediction    = false;
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
    int defaultSubsetSizePerEpoch = 0)
{
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
            }

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
               << std::setw(15) << "MinDelta" << "\n";
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
                   << std::setw(15) << s.curriculumMinDelta[i] << "\n";
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
            p.epsPrediction, initLWP, p.batchSize, initGC, initLR,
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
                s.stage1Model->updateUseDimWeightController(phase0UseDimWeightController);
                mf::LogInfo(moduleName)
                    << "Watch for parameter override: Loaded stage-1 model checkpoint from " << s.ckptStage1File;
            } else {
                s.stage1Model = makeSBDM(3, 0);
            }
            if (s.trainingSize > 0) s.stage1TrainingData.reserve(s.trainingSize);
            else                    s.stage1TrainingData.reserve(1000);
        }
        if (s.trainStage2) {
            if (!s.ckptStage2File.empty()) {
                s.stage2Model = std::make_unique<ScoreBasedDiffusionModel>(
                    ScoreBasedDiffusionModel::loadModel(rf, rg, s.ckptStage2File));
                s.stage2Model->updateUseDimWeightController(phase0UseDimWeightController);
                mf::LogInfo(moduleName)
                    << "Watch for parameter override: Loaded stage-2 model checkpoint from " << s.ckptStage2File;
            } else {
                s.stage2Model = makeSBDM(3, 3);
            }
            if (s.trainingSize > 0) s.stage2TrainingData.reserve(s.trainingSize);
            else                    s.stage2TrainingData.reserve(1000);
        }
    } else {
        if (!s.ckptAllAtOnceFile.empty()) {
            s.allAtOnceModel = std::make_unique<ScoreBasedDiffusionModel>(
                ScoreBasedDiffusionModel::loadModel(rf, rg, s.ckptAllAtOnceFile));
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
inline void accumulateNorm(TrainState& s, double t_trans, double x_trans, double y_trans,
                           double pr_t, double pphi_t, double pz_t)
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
    update(s.pr_mean,   s.pr_M2,   pr_t);
    update(s.pphi_mean, s.pphi_M2, pphi_t);
    update(s.pz_mean,   s.pz_M2,   pz_t);
}

// ---------------------------------------------------------------------------
// collectSample — push a transformed hit into the appropriate training data vector(s)
// ---------------------------------------------------------------------------
inline void collectSample(TrainState& s, double t_trans, double x_trans, double y_trans,
                           double pr_t, double pphi_t, double pz_t)
{
    if (s.useTwoStageTraining) {
        if (s.trainStage1) {
            DiffusionTrainingSample s1;
            s1.x = {t_trans, x_trans, y_trans};
            s.stage1TrainingData.push_back(std::move(s1));
        }
        if (s.trainStage2) {
            DiffusionTrainingSample s2;
            s2.x    = {pr_t, pphi_t, pz_t};
            s2.cond = {t_trans, x_trans, y_trans};
            s.stage2TrainingData.push_back(std::move(s2));
        }
    } else {
        DiffusionTrainingSample sa;
        sa.x = {t_trans, x_trans, y_trans, pr_t, pphi_t, pz_t};
        s.allAtOnceTrainingData.push_back(std::move(sa));
    }
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
//   learned it. Branches: t, then per dimension i: true_dim<i>, denoised_dim<i>,
//   in the model's dimension order (all-at-once: t,x,y,pr,pphi,pz; stage 2:
//   pr,pphi,pz).
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
        tree->Branch("t", &tBranch);
        for (int i = 0; i < dim; ++i) {
            tree->Branch(("true_dim" + std::to_string(i)).c_str(), &trueD[i]);
            tree->Branch(("denoised_dim" + std::to_string(i)).c_str(), &denoisedD[i]);
        }

        size_t written = 0;
        for (size_t k = 0; k < data.size() && written < nUse; k += stride, ++written) {
            const auto& sample = data[k];
            auto res = model.denoiseOneStep(sample.x, sample.cond, t, s.denoiseDiagnosticUseEMA);
            for (int i = 0; i < dim; ++i) {
                trueD[i]     = sample.x[i] * stdev[i] + mean[i];
                denoisedD[i] = res.value[i];
            }
            tree->Fill();
        }
        tree->Write();
    }
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
//   t0 is responsible. Branches: t0, then per dimension i: true_dim<i>,
//   sampled_dim<i>, in the model's dimension order.
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
        tree->Branch("t0", &t0Branch);
        for (int i = 0; i < dim; ++i) {
            tree->Branch(("true_dim" + std::to_string(i)).c_str(), &trueD[i]);
            tree->Branch(("sampled_dim" + std::to_string(i)).c_str(), &sampledD[i]);
        }

        size_t written = 0;
        for (size_t k = 0; k < data.size() && written < nUse; k += stride, ++written) {
            const auto& sample = data[k];
            auto res = model.partialReverseSample(sample.x, sample.cond, t0,
                s.denoiseDiagnosticUseEMA, s.partialReverseUseHeun, s.partialReverseUseSDE,
                s.partialReverseDiffusionSteps, s.partialReverseSdeToOdeSigmaThreshold);
            for (int i = 0; i < dim; ++i) {
                trueD[i]    = sample.x[i] * stdev[i] + mean[i];
                sampledD[i] = res.value[i];
            }
            tree->Fill();
        }
        tree->Write();
        mf::LogInfo(moduleName)
            << "Partial-reverse diagnostic tree " << treeName << " filled (" << written << " samples)";
    }
    fout.Close();
    mf::LogInfo(moduleName)
        << "Partial-reverse diagnostic written to " << outFile
        << " (" << nUse << " samples per t0, " << s.partialReverseT0s.size() << " t0 values)";
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
    s.pr_stdev   = stdev(s.pr_M2);
    s.pphi_stdev = stdev(s.pphi_M2);
    s.pz_stdev   = stdev(s.pz_M2);

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
    const bool diagnosticMode = !s.denoiseDiagnosticTs.empty() || !s.partialReverseT0s.empty();
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
    };

    // Per-epoch training loop shared by all three model paths.
    // resetPhase0: if true, resets currentBias/tLow to phase-0 values before the loop
    //              (needed for stage-2, which follows stage-1 in the same endJob call).
    auto trainLoop = [&](ScoreBasedDiffusionModel& model,
                         std::vector<DiffusionTrainingSample>& data,
                         const std::string& outFile,
                         bool resetPhase0) {
        if (resetPhase0) {
            s.currentBiasLowSigma    = s.nPhase > 1 ? s.curriculumBiasLowSigma[0]    : s.currentBiasLowSigma;
            s.currentTLowBound       = s.nPhase > 1 ? s.curriculumTLowBound[0]       : s.currentTLowBound;
            s.currentTFocusLow       = s.nPhase > 1 ? s.curriculumTFocusLow[0]       : s.currentTFocusLow;
            s.currentTFocusHigh      = s.nPhase > 1 ? s.curriculumTFocusHigh[0]      : s.currentTFocusHigh;
            s.currentTFocusFraction  = s.nPhase > 1 ? s.curriculumTFocusFraction[0]  : s.currentTFocusFraction;
            s.currentSubsetSizePerEpoch = s.nPhase > 1 ? s.curriculumSubsetSizePerEpoch[0] : s.currentSubsetSizePerEpoch;
        }
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
                << ", trainingSubsetSizePerEpoch=" << s.currentSubsetSizePerEpoch;
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
                            s.currentTFocusLow, s.currentTFocusHigh, s.currentTFocusFraction);
                if (std::find(s.saveEpochs.begin(), s.saveEpochs.end(), e) != s.saveEpochs.end()) {
                    model.saveModel(base + ".epoch" + std::to_string(e) + ".bin");
                    if (s.saveAlsoCsv)
                        model.saveModelCsv(base + ".epoch" + std::to_string(e) + ".csv");
                }
            }
            model.saveModel(outFile);
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
                            s.currentTFocusLow, s.currentTFocusHigh, s.currentTFocusFraction);
                double loss = model.getLastEpochLoss(); // read immediately after train()
                const bool conv = tracer.update(loss);

                // Per-epoch planner diagnostics: smoothed loss, its relative improvement
                // vs the running best (the "delta" compared against minDelta), and the
                // no-improvement streak toward convergence. Makes convergence behavior
                // readable from the log instead of reverse-engineered from raw Loss.
                {
                    std::ostringstream pss;
                    pss << "Phase " << phaseNum << " planner: rawLoss=" << loss
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
                    model.saveModel(phaseTag + ".bestRaw.bin");
                    if (s.saveAlsoCsv)
                        model.saveModelCsv(phaseTag + ".bestRaw.csv");
                }
                if (tracer.smoothedImproved) {
                    // Capture the smoothed-best in memory for resume-from-best, and on disk.
                    model.snapshotNetwork();
                    model.saveModel(phaseTag + ".bestSmoothed.bin");
                    if (s.saveAlsoCsv)
                        model.saveModelCsv(phaseTag + ".bestSmoothed.csv");
                }
                // Honor explicit checkpoint epochs (indexed by global epoch count).
                if (std::find(s.saveEpochs.begin(), s.saveEpochs.end(), globalEpoch) != s.saveEpochs.end()) {
                    model.saveModel(base + ".epoch" + std::to_string(globalEpoch) + ".bin");
                    if (s.saveAlsoCsv)
                        model.saveModelCsv(base + ".epoch" + std::to_string(globalEpoch) + ".csv");
                }

                if (conv) { converged = true; break; }
                if (epochInPhase >= maxEpochs) break;
            }

            // Phase-final snapshot: the literal final-epoch state, saved BEFORE any
            // resume-from-best restore so ".final" always means "last epoch trained".
            model.saveModel(phaseTag + ".final.bin");
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
                model.saveModel(outFile);
                if (s.saveAlsoCsv)
                    model.saveModelCsv(base + ".csv");
                return; // abort: surface the bad plan instead of advancing
            }
        }
        model.saveModel(outFile);
        if (s.saveAlsoCsv)
            model.saveModelCsv(base + ".csv");
        mf::LogInfo(moduleName) << "Training completed (auto planner), saved to " << outFile;
    };

    if (s.useTwoStageTraining) {
        if (s.trainStage1) {
            if (s.trainingSize > 0 && (int)s.stage1TrainingData.size() > s.trainingSize)
                s.stage1TrainingData.resize(s.trainingSize);
            s.stage1Model->normalizeData(
                {s.t_mean, s.x_mean, s.y_mean},
                {s.t_stdev, s.x_stdev, s.y_stdev},
                s.stage1TrainingData);
            if (diagnosticMode) {
                runDiagnostics(*s.stage1Model, s.stage1TrainingData,
                    {s.t_mean, s.x_mean, s.y_mean},
                    {s.t_stdev, s.x_stdev, s.y_stdev},
                    s.stage1ModelFile);
            } else {
                mf::LogInfo(moduleName)
                    << "Training stage-1 diffusion model with " << s.stage1TrainingData.size()
                    << " samples and " << s.trainingEpochs << " epochs...";
                trainLoop(*s.stage1Model, s.stage1TrainingData, s.stage1ModelFile, false);
            }
        }
        if (s.trainStage2) {
            if (s.trainingSize > 0 && (int)s.stage2TrainingData.size() > s.trainingSize)
                s.stage2TrainingData.resize(s.trainingSize);
            s.stage2Model->normalizeData(
                {s.pr_mean, s.pphi_mean, s.pz_mean, s.t_mean, s.x_mean, s.y_mean},
                {s.pr_stdev, s.pphi_stdev, s.pz_stdev, s.t_stdev, s.x_stdev, s.y_stdev},
                s.stage2TrainingData);
            if (diagnosticMode) {
                runDiagnostics(*s.stage2Model, s.stage2TrainingData,
                    {s.pr_mean, s.pphi_mean, s.pz_mean, s.t_mean, s.x_mean, s.y_mean},
                    {s.pr_stdev, s.pphi_stdev, s.pz_stdev, s.t_stdev, s.x_stdev, s.y_stdev},
                    s.stage2ModelFile);
            } else {
                mf::LogInfo(moduleName)
                    << "Training stage-2 diffusion model with " << s.stage2TrainingData.size()
                    << " samples and " << s.trainingEpochs << " epochs...";
                trainLoop(*s.stage2Model, s.stage2TrainingData, s.stage2ModelFile, true);
            }
        }
    } else {
        if (s.allAtOnceModelFile.empty())
            throw cet::exception(moduleName) << "All-at-once training requires SBDMallAtOnceModelFile.";
        if (s.trainingSize > 0 && (int)s.allAtOnceTrainingData.size() > s.trainingSize)
            s.allAtOnceTrainingData.resize(s.trainingSize);
        s.allAtOnceModel->normalizeData(
            {s.t_mean, s.x_mean, s.y_mean, s.pr_mean, s.pphi_mean, s.pz_mean},
            {s.t_stdev, s.x_stdev, s.y_stdev, s.pr_stdev, s.pphi_stdev, s.pz_stdev},
            s.allAtOnceTrainingData);
        if (diagnosticMode) {
            runDiagnostics(*s.allAtOnceModel, s.allAtOnceTrainingData,
                {s.t_mean, s.x_mean, s.y_mean, s.pr_mean, s.pphi_mean, s.pz_mean},
                {s.t_stdev, s.x_stdev, s.y_stdev, s.pr_stdev, s.pphi_stdev, s.pz_stdev},
                s.allAtOnceModelFile);
        } else {
            mf::LogInfo(moduleName)
                << "Training all-at-once diffusion model with " << s.allAtOnceTrainingData.size()
                << " samples and " << s.trainingEpochs << " epochs...";
            trainLoop(*s.allAtOnceModel, s.allAtOnceTrainingData, s.allAtOnceModelFile, false);
        }
    }
}

} // namespace VDResampler
} // namespace mu2e
