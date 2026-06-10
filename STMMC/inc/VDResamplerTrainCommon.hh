#pragma once

// Shared state and logic for VDResamplerTrain and VDResamplerTrainFromRoot modules.
// Header-only — all functions are inline. Include pattern follows VDResamplerTransforms.hh.
// Yongyi Wu, Jun. 2026

#include <algorithm>
#include <cmath>
#include <iomanip>
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
    int trainingSubsetSizePerEpoch  = 0;
    std::vector<int> saveEpochs;
    bool saveAlsoCsv                = false; // if true, also write a CSV alongside every binary save

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
    bool   promoteEMAOnStart   = false;
    bool   currentBiasLowSigma = false;
    double currentTLowBound    = 0.0;
    std::vector<int> phaseBoundaries;

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
// ModelBuildParams — SBDM constructor arguments, filled from fhicl Config
// ---------------------------------------------------------------------------
struct ModelBuildParams {
    int    timeEmbeddingDim = 0;
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
    bool defaultPromoteEMA = false)
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
               << std::setw(15) << "PromoteEMA" << "\n";
            for (int i = 0; i < s.nPhase; ++i)
                ss << std::setw(10) << s.curriculumEpochs[i]
                   << std::setw(20) << s.curriculumLossWeightPower[i]
                   << std::setw(20) << s.curriculumGradientClip[i]
                   << std::setw(20) << s.curriculumLearningRate[i]
                   << std::setw(15) << s.curriculumBiasLowSigma[i]
                   << std::setw(15) << s.curriculumTLowBound[i]
                   << std::setw(15) << s.curriculumBatchSize[i]
                   << std::setw(15) << s.curriculumPromoteEMA[i] << "\n";
            ss << "[End of Curriculum Training Schema]";
            mf::LogInfo(moduleName) << ss.str();
        }
    }
    s.currentBiasLowSigma = s.nPhase > 1 ? s.curriculumBiasLowSigma[0] : defaultBiasLowSigma;
    s.currentTLowBound    = s.nPhase > 1 ? s.curriculumTLowBound[0]    : defaultTLowBound;
}

// ---------------------------------------------------------------------------
// buildModels — construct SBDM model objects (or load checkpoints) and reserve
//               training data vectors.  Sets trainStage1/trainStage2 flags.
// ---------------------------------------------------------------------------
inline void buildModels(TrainState& s, const ModelBuildParams& p,
                        CLHEP::RandFlat& rf, CLHEP::RandGaussQ& rg,
                        const std::string& moduleName)
{
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
        return std::make_unique<ScoreBasedDiffusionModel>(
            rf, rg, dim, condDim, p.timeEmbeddingDim, p.hidden, p.layers, p.optimizer,
            p.adamBeta1, p.adamBeta2, p.adamEps, p.noiseSchedule,
            p.betaMin, p.betaMax, p.cosineOffset, p.logSigMin, p.logSigMax,
            p.epsPrediction, initLWP, p.batchSize, initGC, initLR,
            p.useDimWeightController, p.dimWeightControllerEMADecay,
            p.useEMANetwork, p.emaNetworkDecay, p.diffusionSteps
        );
    };

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
// runTraining — finalise normalisation, then run the full per-epoch training
//               loop (with curriculum phase transitions and epoch checkpointing)
//               for each active model path.
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

    // Per-epoch training loop shared by all three model paths.
    // resetPhase0: if true, resets currentBias/tLow to phase-0 values before the loop
    //              (needed for stage-2, which follows stage-1 in the same endJob call).
    auto trainLoop = [&](ScoreBasedDiffusionModel& model,
                         std::vector<DiffusionTrainingSample>& data,
                         const std::string& outFile,
                         bool resetPhase0) {
        if (resetPhase0) {
            s.currentBiasLowSigma = s.nPhase > 1 ? s.curriculumBiasLowSigma[0] : s.currentBiasLowSigma;
            s.currentTLowBound    = s.nPhase > 1 ? s.curriculumTLowBound[0]    : s.currentTLowBound;
        }
        if (s.promoteEMAOnStart || (s.nPhase > 1 && s.curriculumPromoteEMA[0]))
            model.promoteEMAToNetwork();
        mf::LogInfo(moduleName)
            << "Initial training settings: biasLowSigma=" << s.currentBiasLowSigma
            << ", tLowBound=" << s.currentTLowBound
            << ", trainingSubsetSizePerEpoch=" << s.trainingSubsetSizePerEpoch;
        const std::string base = stripExt(outFile);
        for (int e = 1; e <= s.trainingEpochs; ++e) {
            if (s.nPhase > 1) {
                for (int ph = 0; ph < s.nPhase - 1; ++ph) {
                    if (e == s.phaseBoundaries[ph] + 1) {
                        double lwp = s.curriculumLossWeightPower[ph + 1];
                        double gc  = s.curriculumGradientClip   [ph + 1];
                        double lr  = s.curriculumLearningRate   [ph + 1];
                        int    bs  = s.curriculumBatchSize      [ph + 1];
                        s.currentBiasLowSigma = s.curriculumBiasLowSigma[ph + 1];
                        s.currentTLowBound    = s.curriculumTLowBound   [ph + 1];
                        mf::LogInfo(moduleName)
                            << "Switching to phase " << (ph + 2)
                            << ": lossWeightPower=" << lwp
                            << ", gradientClipThreshold=" << gc
                            << ", learningRate=" << lr
                            << ", batchSize=" << bs
                            << ", biasLowSigma=" << s.currentBiasLowSigma
                            << ", tLowBound=" << s.currentTLowBound
                            << ", trainingSubsetSizePerEpoch=" << s.trainingSubsetSizePerEpoch;
                        model.updateLossWeightPower(lwp);
                        model.updateGradientClipThreshold(gc);
                        model.updateLearningRate(lr);
                        model.updateBatchSize(bs);
                        if (s.curriculumPromoteEMA[ph + 1])
                            model.promoteEMAToNetwork();
                    }
                }
            }
            mf::LogInfo(moduleName) << "Epoch " << e << "/" << s.trainingEpochs;
            model.train(data, 1, s.trainingSubsetSizePerEpoch, s.currentBiasLowSigma, s.currentTLowBound);
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
    };

    if (s.useTwoStageTraining) {
        if (s.trainStage1) {
            if (s.trainingSize > 0 && (int)s.stage1TrainingData.size() > s.trainingSize)
                s.stage1TrainingData.resize(s.trainingSize);
            s.stage1Model->normalizeData(
                {s.t_mean, s.x_mean, s.y_mean},
                {s.t_stdev, s.x_stdev, s.y_stdev},
                s.stage1TrainingData);
            mf::LogInfo(moduleName)
                << "Training stage-1 diffusion model with " << s.stage1TrainingData.size()
                << " samples and " << s.trainingEpochs << " epochs...";
            trainLoop(*s.stage1Model, s.stage1TrainingData, s.stage1ModelFile, false);
        }
        if (s.trainStage2) {
            if (s.trainingSize > 0 && (int)s.stage2TrainingData.size() > s.trainingSize)
                s.stage2TrainingData.resize(s.trainingSize);
            s.stage2Model->normalizeData(
                {s.pr_mean, s.pphi_mean, s.pz_mean, s.t_mean, s.x_mean, s.y_mean},
                {s.pr_stdev, s.pphi_stdev, s.pz_stdev, s.t_stdev, s.x_stdev, s.y_stdev},
                s.stage2TrainingData);
            mf::LogInfo(moduleName)
                << "Training stage-2 diffusion model with " << s.stage2TrainingData.size()
                << " samples and " << s.trainingEpochs << " epochs...";
            trainLoop(*s.stage2Model, s.stage2TrainingData, s.stage2ModelFile, true);
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
        mf::LogInfo(moduleName)
            << "Training all-at-once diffusion model with " << s.allAtOnceTrainingData.size()
            << " samples and " << s.trainingEpochs << " epochs...";
        trainLoop(*s.allAtOnceModel, s.allAtOnceTrainingData, s.allAtOnceModelFile, false);
    }
}

} // namespace VDResampler
} // namespace mu2e
