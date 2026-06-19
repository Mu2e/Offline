// Module for training and using (variance preserving) score-based diffusion model
// Added by Yongyi Wu
// Mar. 2026

#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>
#include <cassert>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

namespace mu2e{
    struct DiffusionTrainingSample {
        std::vector<double> x;    // transformed state vector to diffuse (size = dim)
        std::vector<double> cond; // optional conditioning vector (size = conditionDim)
    };

    struct SBDMGeneratedSample {
        std::vector<double> zscore; // normalized data
        std::vector<double> value;  // unnormalized data
    };

    // One peak-importance-sampling window (see ScoreBasedDiffusionModel::train). Oversamples
    // training examples whose NORMALIZED coordinate `dim` lies in [low, high) to fix the
    // under-weighting of a rare but important localized feature. Multiple windows may be
    // supplied; they should be disjoint (a sample is assigned to the first window it matches).
    struct PeakWindow {
        int    dim    = -1;   // normalized state-coordinate index this window is defined on
        double low    = 0.0;  // window lower edge (z-score units), inclusive
        double high   = 0.0;  // window upper edge (z-score units), exclusive
        double gMax   = 0.0;  // peak sampling-fraction ceiling at low sigma (0 < gMax < 1; sum over windows < 1)
        double sigma0 = 0.0;  // Gaussian sigma-taper scale (~ feature width); oversampling decays for sigma >> sigma0
        double alpha  = 1.0;  // 1 = unbiased (variance reduction only), <1 = deliberate up-weight
    };

    // One conditional-loss diagnostic sample: the per-dimension epsilon-prediction squared error
    // for one data point at a freshly drawn diffusion time t. Used to profile L(sigma) inside vs
    // outside a feature window (see ScoreBasedDiffusionModel::evalEpsLossSample).
    struct SBDMEpsLossSample {
        double t = 0.0;                 // drawn diffusion time
        double sigma = 0.0;             // sigma(t)
        std::vector<double> perDimLoss; // (eps_hat - eps)^2 per state dimension (size dim)
    };

    // First-layer input-feature-block magnitude (see ScoreBasedDiffusionModel::firstLayerBlockMagnitudes).
    // Reports, for one contiguous block of the network input, the L2 norm of the first layer's
    // weight columns and of the (training-gradient) gradient columns over that block.
    struct SBDMFeatureBlockMagnitude {
        std::string name;     // human-readable block label, e.g. "fourier_state[5]"
        int    kind  = 0;     // 0 raw-state, 1 fourier-state, 2 raw-cond, 3 fourier-cond, 4 raw-time, 5 fourier-time
        int    coord = -1;    // state/condition coordinate index this block belongs to (-1 for time)
        int    nCols = 0;     // number of input columns in the block
        double weightL2 = 0.0;
        double gradL2   = 0.0;
    };

    class ScoreBasedDiffusionModel {
    public:
        // Enumeration for optimizer selection
        enum class OptimizerType {
            SGD,   // Stochastic Gradient Descent
            ADAM   // Adam optimizer
        };

        // Enumeration for noise schedule selection
        enum class NoiseScheduleType {
            LINEAR,  // Linear noise schedule
            COSINE,  // Cosine noise schedule
            LOGSIG   // Log-sigma (pure exponential) noise schedule: sigma(t) = sigMin*exp(k*t), k=ln(sigMax/sigMin)
        };

        // Constructor: Initialize diffusion model with CLHEP random distributions.
        //
        // Parameters:
        //   randFlat                - Reference to CLHEP RandFlat for uniform sampling (externally managed)
        //   randGaussQ              - Reference to CLHEP RandGaussQ for Gaussian noise (externally managed)
        //   dim                     - Dimensionality of the state space
        //   conditionDim            - Dimensionality of the optional conditioning vector (default: 0 for unconditional model)
        //   timeEmbeddingDim        - Dimensionality of the sinusoidal time embedding applied to the diffusion time t in [0,1].
        //                             0 = raw scalar t (default, backward-compatible).
        //                             Even integer >= 2 = sinusoidal embedding: pairs [sin(2π·2^i·t), cos(2π·2^i·t)] for i=0..k/2-1.
        //   inputEmbeddingDims      - Per-coordinate depth of the sinusoidal (Fourier) embedding applied to the dim state
        //                             coordinates. One entry per state dimension. Accepts {} (no embedding on any dim),
        //                             {k} (broadcast depth k to every dim), or a length-dim vector (per-dim depth).
        //                             Each depth is 0 (raw coordinate only) or an even integer >= 2, emitting pairs
        //                             [sin(π·2^i·x_j), cos(π·2^i·x_j)] for i=0..k/2-1 on that coordinate.
        //                             Counteracts MLP spectral bias so the score network can represent structure much finer
        //                             than O(1) in the normalized coordinates (e.g. a narrow peak). To resolve a feature of
        //                             normalized width w, the top frequency π·2^(k/2-1) should be ≳ 2π/w (k=24 reaches w~1e-3).
        //                             Use per-dim depth to give a structured coordinate (e.g. pz) a deep embedding while
        //                             smooth coordinates get 0, avoiding injecting noise features on dims that don't need them.
        //   conditionEmbeddingDims  - Same per-coordinate Fourier embedding, applied to the conditionDim condition
        //                             coordinates ({} / {k} / length-conditionDim). 0 on every dim unless the conditional
        //                             distribution changes sharply with that condition coordinate. Must be empty when conditionDim==0.
        //   hidden                  - Size of hidden layers in the neural network
        //   layers                  - Number of layers in the network
        //   optimizerType           - Type of optimizer to use (SGD or ADAM, default: ADAM)
        //   adamBeta1               - Adam optimizer beta1 parameter (default: 0.9)
        //   adamBeta2               - Adam optimizer beta2 parameter (default: 0.999)
        //   adamEps                 - Adam optimizer epsilon parameter (default: 1e-8)
        //   scheduleType            - Type of noise schedule (LINEAR, COSINE, or LOGSIG, default: COSINE)
        //   betaMin                 - Minimum noise schedule parameter (for LINEAR schedule, default: 1e-4)
        //   betaMax                 - Maximum noise schedule parameter (for LINEAR schedule, default: 0.02)
        //   cosineOffset            - Offset parameter (for cosine schedule, default: 0.008)
        //   logSigMin               - Minimum sigma for LOGSIG schedule (default: 1e-5)
        //   logSigMax               - Maximum sigma for LOGSIG schedule (default: 1.0)
        //   epsPrediction           - If true, network predicts noise epsilon instead of the score (default: false)
        //   lossWeightPower         - Power for weighting the loss function (default: 2.0 for quadratic weighting)
        //   batchSize               - Batch size for training (default: 32)
        //   gradientClipThreshold   - Threshold for gradient clipping (default: 1.0)
        //   learningRate            - Learning rate for training (default: 1e-3)
        //   useDimWeightController  - If true, enables adaptive per-dimension gradient weighting (default: false)
        //   dimWeightEMADecay       - EMA decay rate for per-dimension loss tracking (default: 0.99)
        //   useEMANetwork           - If true, maintains a slow-moving EMA copy of network weights for inference (default: true)
        //   emaNetworkDecay         - Decay rate for EMA network update per optimizer step (default: 0.9999)
        //   diffusionSteps          - Number of steps in the diffusion process (default: 200)
        //   initializeRandomWeights - If true, initialize network weights from random Gaussian draws.
        //                             Set false when constructing (loading) from a saved model.
        ScoreBasedDiffusionModel(
            // Network architecture parameters
            CLHEP::RandFlat& randFlat,
            CLHEP::RandGaussQ& randGaussQ,
            int dim,
            int conditionDim,
            int timeEmbeddingDim = 0,
            std::vector<int> inputEmbeddingDims = {},
            std::vector<int> conditionEmbeddingDims = {},
            int hidden = 128,
            int layers = 4,
            // Optimizer configuration
            OptimizerType optimizerType = OptimizerType::ADAM,
            double adamBeta1 = 0.9,
            double adamBeta2 = 0.999,
            double adamEps = 1e-8,
            // Noise schedule configuration
            NoiseScheduleType scheduleType = NoiseScheduleType::COSINE,
            // -- linear schedule parameters
            double betaMin = 1e-4,
            double betaMax = 0.02,
            // -- cosine schedule parameters
            double cosineOffset = 0.008,
            // -- log-sigma schedule parameters
            double logSigMin = 1e-5,
            double logSigMax = 1.0,
            // Training configuration
            // -- Training target
            bool epsPrediction = false,
            // -- Training configuration
            double lossWeightPower = 2.0,
            int batchSize = 32,
            double gradientClipThreshold = 1.0,
            double learningRate = 1e-3,
            // -- Adaptive dimensional weight controller
            bool useDimWeightController = false,
            double dimWeightEMADecay = 0.99,
            // EMA copy of network parameters for inference
            bool useEMANetwork = true,
            double emaNetworkDecay = 0.9999,
            // Diffusion process configuration
            int diffusionSteps = 200,
            bool initializeRandomWeights = true
        );

        // Data normalization to be applied before training
        // Parameters:
        //   mean - Mean values for each dimension of the data (both state and condition)
        //   stdev - Standard deviation values for each dimension of the data (both state and condition)
        //   data - Training samples (transformed state vectors)
        void normalizeData(
            const std::vector<double>& mean,
            const std::vector<double>& stdev,
            std::vector<DiffusionTrainingSample>& data
        );

        // functions to update certain training parameters
        double updateLossWeightPower(
            double value
        ){
            lossWeightPower_ = value;
            return lossWeightPower_;
        }

        double updateGradientClipThreshold(
            double value
        ){
            gradientClipThreshold_ = value;
            return gradientClipThreshold_;
        }

        double updateLearningRate(
            double value
        ){
            learningRate_ = value;
            return learningRate_;
        }

        int updateBatchSize(
            int value
        ){
            batchSize_ = value;
            if (useEMANetwork_) {
                emaNetworkDecay_ = std::pow(emaNetworkDecayBase_,
                                           (double)batchSize_ / kEMABatchSizeRef_);
                mf::LogInfo("ScoreBasedDiffusionModel")
                    << "Batch size updated to " << batchSize_
                    << "; EMA decay rescaled to " << emaNetworkDecay_;
            }
            return batchSize_;
        }

        // Enable/disable the adaptive per-dimension gradient weight controller mid-run
        // (e.g. as a curriculum-phase change, or as an explicit override after loading a
        // checkpoint that was trained with it on). Turning it OFF does NOT reset
        // dimWeights_: train() keeps multiplying each dimension's gradient by the frozen
        // weight, it merely stops adapting them. The frozen values are logged so the
        // override is auditable. Returns the resulting flag state.
        bool updateUseDimWeightController(
            bool enabled
        ){
            if (useDimWeightController_ && !enabled) {
                std::ostringstream woss;
                woss << "Dimensional weight controller turned OFF; freezing dimWeights at [";
                for (int i = 0; i < dim_; ++i) {
                    woss << dimWeights_[i];
                    if (i < dim_ - 1) woss << ", ";
                }
                woss << "] (still applied to gradients, but no longer adapting).";
                mf::LogInfo("ScoreBasedDiffusionModel::updateUseDimWeightController") << woss.str();
            }
            useDimWeightController_ = enabled;
            return useDimWeightController_;
        }

        void promoteEMAToNetwork() {
            if (!useEMANetwork_) {
                mf::LogWarning("ScoreBasedDiffusionModel")
                    << "promoteEMAToNetwork called but EMA network is disabled — no-op.";
                return;
            }
            for (size_t l = 0; l < network_.size(); ++l) {
                network_[l].W = emaNetwork_[l].W;
                network_[l].b = emaNetwork_[l].b;
                for (auto& row : network_[l].mW) std::fill(row.begin(), row.end(), 0.0);
                for (auto& row : network_[l].vW) std::fill(row.begin(), row.end(), 0.0);
                std::fill(network_[l].mb.begin(), network_[l].mb.end(), 0.0);
                std::fill(network_[l].vb.begin(), network_[l].vb.end(), 0.0);
            }
            adamStep_ = 0;
            mf::LogInfo("ScoreBasedDiffusionModel::promoteEMAToNetwork")
                << "EMA weights promoted to network. Adam optimizer state reset.";
        }

        const std::vector<double>& getDimWeights() const { return dimWeights_; }

        // In-memory snapshot of the trainable state. Used by the auto curriculum planner to
        // capture the smoothed-best point within a phase and restore it before advancing —
        // the substitute for EMA promotion when that is disabled. The snapshot captures the
        // SAME mutable state that loadModel() round-trips, so an in-memory restore is
        // equivalent to reloading the checkpoint: network weights + Adam moments (full
        // Layer), the EMA copy, adamStep_, AND the dimension-weight controller state
        // (dimWeights_ / dimLossEMA_). Omitting the controller state previously left a
        // skewed dimWeights_ in force after a restore, which destabilised the next phase.
        // The model cannot be reassigned in place (reference members delete operator=), so a
        // value-copy snapshot is used instead of a loadModel() round-trip.
        void snapshotNetwork() {
            networkSnapshot_    = network_;
            emaNetworkSnapshot_ = emaNetwork_;
            adamStepSnapshot_   = adamStep_;
            dimWeightsSnapshot_ = dimWeights_;
            dimLossEMASnapshot_ = dimLossEMA_;
            hasSnapshot_        = true;
        }
        void restoreNetwork() {
            if (!hasSnapshot_) {
                mf::LogWarning("ScoreBasedDiffusionModel")
                    << "restoreNetwork called with no snapshot — no-op.";
                return;
            }
            network_    = networkSnapshot_;
            emaNetwork_ = emaNetworkSnapshot_;
            adamStep_   = adamStepSnapshot_;
            dimWeights_ = dimWeightsSnapshot_;
            dimLossEMA_ = dimLossEMASnapshot_;
            mf::LogInfo("ScoreBasedDiffusionModel::restoreNetwork")
                << "Restored network weights, optimizer, and dim-weight controller from in-memory snapshot.";
        }
        bool hasNetworkSnapshot() const { return hasSnapshot_; }
        void clearNetworkSnapshot() {
            networkSnapshot_.clear();
            emaNetworkSnapshot_.clear();
            dimWeightsSnapshot_.clear();
            dimLossEMASnapshot_.clear();
            hasSnapshot_ = false;
        }

        // Train the score network on a batch of samples.
        // Uses random sampling and noise injection via the external engine.
        // Note that training needs to occur on all data samples. Training on multiple small subsets
        // of data and then averaging or aggregating the model parameters is not supported and may lead to
        // issues as neural networks are not linear.
        //
        // Parameters:
        //   data         - Training samples (transformed state vectors)
        //   epochs       - Number of training epochs to perform
        //   trainSubsetDataSize - If > 0, randomly samples this many data points from the full dataset for each epoch (default: 0, uses all data)
        //   biasLowSigma  - If true, biases the sampling of diffusion time t towards smaller values (smaller sigma) by sampling t^2 instead of t (default: false)
        //   tLowBound     - If > 0, enforces a lower bound on the sampled diffusion time t to focus training on larger sigma values (default: 0.0, no lower bound)
        //   tFocusLow/tFocusHigh/tFocusFraction
        //                 - If tFocusFraction > 0, each training sample draws t uniformly from the window
        //                   [tFocusLow, tFocusHigh] with probability tFocusFraction, and from the regular
        //                   tLowBound/biasLowSigma logic otherwise. Concentrates gradient steps on a target
        //                   sigma band (e.g. the band matching a narrow spectral feature) while keeping full
        //                   [0,1] coverage. The t-sampling distribution only reweights the loss; it does not
        //                   bias the learned data distribution. (defaults: 0.0, 0.0, 0.0 = disabled)
        //   peakWindows   - Peak importance sampling (disabled when empty). A list of disjoint windows
        //                   (see PeakWindow); each oversamples training examples whose NORMALIZED
        //                   coordinate `dim` lies in [low, high) to fix the under-weighting of rare but
        //                   physically important features (e.g. narrow pz peaks holding a tiny data
        //                   fraction). Per draw, t is sampled first; window k is drawn with probability
        //                   g_eff_k(t) = max(f_k, gMax_k * exp(-sigma(t)^2 / (2 sigma0_k^2))), where f_k is
        //                   the empirical in-window fraction — so oversampling concentrates at
        //                   sigma <~ sigma0_k (set ~ the feature width) and decays to the natural rate at
        //                   high sigma. A single cumulative draw selects window k or the out-of-window
        //                   pool. The loss is reweighted by (f_k/g_eff_k)^alpha_k in window k (alpha_k = 1
        //                   exactly unbiased / variance reduction only; alpha_k < 1 a deliberate up-weight),
        //                   and by (f_Q/(1-Sum g_eff)) out (unbiased continuum), with f_Q = 1 - Sum f_k.
        //                   Requires Sum gMax_k < 1 (leaves probability for the out-of-window pool). Each
        //                   pool is drawn without replacement via a per-epoch reshuffled cursor.
        void train(
            const std::vector<DiffusionTrainingSample>& data,
            int epochs,
            int trainSubsetDataSize = 0,
            bool biasLowSigma = false,
            double tLowBound = 0.0,
            double tFocusLow = 0.0,
            double tFocusHigh = 0.0,
            double tFocusFraction = 0.0,
            const std::vector<PeakWindow>& peakWindows = {}
        );

        // Diagnostic: average per-sample loss at the CURRENT weights with NO optimizer step.
        // Replicates train()'s per-sample loss (same t sampling / focus / bias, addNoise,
        // forward, computeLoss with sigma^lossWeightPower weighting) over a subset, but does
        // NOT backprop, step the optimizer, update the EMA network, touch dimLossEMA_, or
        // increment adamStep_. Used to read the loss of weights as just installed by
        // restoreNetwork()/promoteEMAToNetwork(), isolating weight quality from the
        // destabilising effect of the first training step. Draws from the RNG (t + noise).
        double evaluateAverageLoss(
            const std::vector<DiffusionTrainingSample>& data,
            int subsetSize = 0,
            bool biasLowSigma = false,
            double tLowBound = 0.0,
            double tFocusLow = 0.0,
            double tFocusHigh = 0.0,
            double tFocusFraction = 0.0
        );

        // Public accessor for the noise level sigma(t) (pure function of the schedule), so
        // diagnostics can set a log-sigma axis without duplicating the schedule.
        double diffusionSigma(double t) const { return sigma(t); }

        // Diagnostic (conditional-loss profile): draw a uniform diffusion time t and fresh noise,
        // run one forward pass, and return the per-dimension epsilon-prediction squared error
        // (eps_hat - eps)^2 along with t and sigma(t). eps_hat is recovered from the network output
        // regardless of epsPrediction_ mode, so the error is comparable across sigma. Binning these
        // by log sigma, split by whether the sample lies inside a feature window, gives L_P(sigma)
        // vs L_Q(sigma): an underfit only at low sigma points to sampling/under-weighting, an
        // underfit at all sigma to model capacity / Fourier resolution. Draws from the RNG; no step.
        SBDMEpsLossSample evalEpsLossSample(
            const std::vector<double>& xNorm,
            const std::vector<double>& condition,
            bool useEMANetworkIfAvailable = true);

        // Diagnostic (first-layer feature-block magnitudes): for each contiguous block of the
        // network input (each raw state coord, each state coord's Fourier columns, each raw
        // condition coord, each condition coord's Fourier columns, and the time block), return the
        // L2 norm of the first layer's weight columns and of the gradient accumulated over nSamples
        // training-like draws (uniform t, addNoise, the configured weighted eps/score loss,
        // backward) with NO optimizer step. Near-zero gradL2 (or weightL2 stuck at init scale) on a
        // block means that input feature is not being used/driven — e.g. dead pz Fourier columns.
        // perDimLossOut is filled with the FRESH mean per-output-dimension squared residual measured
        // over the same sweep (a live alternative to the checkpoint's stale dimLossEMA_). Operates
        // on the base network_ and leaves the gradient buffers zeroed.
        std::vector<SBDMFeatureBlockMagnitude> firstLayerBlockMagnitudes(
            const std::vector<DiffusionTrainingSample>& data,
            int nSamples,
            std::vector<double>& perDimLossOut);

        // Generate a new sample from the diffusion model via reverse process.
        // Uses the external random engine for noise generation during sampling.
        //
        // Parameters:
        //   condition       - Optional conditioning vector (must match conditionDim_ when enabled)
        //   useEMANetworkIfAvailable - If true (default), uses the EMA network when available (i.e. when the model was
        //                     configured with useEMANetwork=true). Pass false to force the base score network.
        //   useHeun         - If true, uses Heun's method (2nd order, default). If false, uses Euler's method (1st order)
        //   useSDE          - If true, uses SDE (Stochastic Differential Equation) method. If false, uses the deterministic reverse process
        //   diffusionSteps  - Number of diffusion steps for sampling (default: -1 uses the model's configured diffusionSteps_)
        //
        // Returns: Two generated sample vectors, zscore and value of dimensions dim_
        SBDMGeneratedSample generateSample(
            const std::vector<double>& condition = {},
            bool useEMANetworkIfAvailable = true,
            bool useHeun = true,
            bool useSDE = true,
            int diffusionSteps = -1,
            double sdeToOdeSigmaThreshold = -1.0
        );

        // One-step denoising diagnostic: perturb a normalized data sample at a fixed diffusion
        // time t, run a single network forward pass, and reconstruct the denoised estimate
        //   x0_hat = (x_t - sigma(t) * eps_hat) / sqrt(alphabar(t)).
        // Comparing the x0_hat distribution against truth separates "the score network never
        // learned a feature" (absent here) from "the reverse-process sampler destroys it"
        // (present here but absent in generated samples).
        //
        // Parameters:
        //   xNorm     - Normalized (z-scored) state vector of size dim_, e.g. a training sample
        //               after normalizeData()
        //   condition - Normalized conditioning vector (must match conditionDim_)
        //   t         - Fixed diffusion time in (0,1); choose sigma(t) comparable to the feature
        //               width under study
        //   useEMANetworkIfAvailable - If true (default), uses the EMA network when available
        //
        // Returns: zscore = x0_hat (normalized space), value = de-normalized x0_hat
        SBDMGeneratedSample denoiseOneStep(
            const std::vector<double>& xNorm,
            const std::vector<double>& condition,
            double t,
            bool useEMANetworkIfAvailable = true
        );

        // Partial-reverse diagnostic: perturb a normalized data sample at diffusion time t0,
        // then run the full multi-step reverse sampler from t0 down to 0 (instead of starting
        // from pure noise at t=1). Scanning t0 localizes where the sampler loses a feature:
        // if the feature survives a start at t0 but not a full generation, the t>t0 phase
        // delivers a wrong marginal; if it is already lost starting at t0, the score or its
        // integration below t0 is responsible. t0=1.0 degenerates to full generation.
        //
        // Parameters:
        //   xNorm     - Normalized (z-scored) state vector of size dim_
        //   condition - Normalized conditioning vector (must match conditionDim_)
        //   t0        - Diffusion time in (0,1] to noise the sample to; snapped to the
        //               sampler's discrete time grid (round(t0*steps)/steps)
        //   useEMANetworkIfAvailable / useHeun / useSDE / diffusionSteps /
        //   sdeToOdeSigmaThreshold - identical to generateSample()
        //
        // Returns: zscore = sampled state (normalized space), value = de-normalized state
        SBDMGeneratedSample partialReverseSample(
            const std::vector<double>& xNorm,
            const std::vector<double>& condition,
            double t0,
            bool useEMANetworkIfAvailable = true,
            bool useHeun = true,
            bool useSDE = true,
            int diffusionSteps = -1,
            double sdeToOdeSigmaThreshold = -1.0
        );

        // Save the model to a binary file (.bin) preserving full double precision.
        // This is the default save format. Use saveModelCsv for human-readable output.
        //
        // Parameters:
        //   filename - Path to the binary file (default: "DiffusionModel.bin")
        void saveModel(const std::string& filename = "DiffusionModel.bin");

        // Save the model parameters to a CSV file with annotations for human inspection.
        //
        // Parameters:
        //   filename - Path to the CSV file (default: "DiffusionModel.csv")
        void saveModelCsv(const std::string& filename = "DiffusionModel.csv");

        // Load model parameters from a file to restore a previously trained model.
        // Auto-detects format by extension: ".bin" loads binary, anything else loads CSV.
        // If optimizer state is present, Adam moments and step counter are restored so
        // training can resume seamlessly.
        //
        // Parameters:
        //   randFlat / RandGaussQ - CLHEP random number generator wrappers being passed
        //   filename              - Path to the file from which model parameters will be loaded
        static ScoreBasedDiffusionModel loadModel(
            CLHEP::RandFlat& randFlat,
            CLHEP::RandGaussQ& randGaussQ,
            const std::string& filename
        );

        // Most recent per-epoch average loss (last entry appended by train()); NaN if
        // epochLosses_ is empty. Each train(..., epochs=1, ...) call appends exactly one
        // entry, so this returns that epoch's loss WHEN CALLED IMMEDIATELY AFTER train().
        // Note: a checkpoint loaded via loadModel() pre-fills epochLosses_ with its saved
        // history, so back() is meaningful as "the current run's latest epoch" only when
        // read right after a fresh train() call (which is how the curriculum planner in
        // VDResamplerTrainCommon uses it).
        double getLastEpochLoss() const {
            return epochLosses_.empty() ? std::numeric_limits<double>::quiet_NaN()
                                        : epochLosses_.back();
        }

    private:

        // ----- network -----

        // Represents a single fully-connected layer with weights and biases.
        struct Layer {
            std::vector<std::vector<double>> W; // Weight matrix [output_size][input_size]
            std::vector<double> b;              // Bias vector [output_size]

            // storage for gradients during back-propagation
            std::vector<std::vector<double>> gradW; // gradient of loss w.r.t. weights
            std::vector<double> gradb;              // gradient of loss w.r.t. biases

            // Adam optimizer state buffers
            std::vector<std::vector<double>> mW; // First moment estimates for weights, m_t = beta1 m_{t-1} + (1-beta1) g_t
            std::vector<std::vector<double>> vW; // Second moment estimates for weights, v_t = beta2 v_{t-1} + (1-beta2) g_t^2

            std::vector<double> mb; // First moment estimates for biases
            std::vector<double> vb; // Second moment estimates for biases
        };

        std::vector<Layer> network_; // Network layers in forward order

        // Storage for activations and pre-activations during forward pass (used in back-propagation)
        // During forward pass, preactivations_[l] = W[l] * activations_[l-1] + b[l],
        //                     and activations_[l] = activationFunction(preactivations_[l])
        // These are needed to compute gradients during the backward pass.
        std::vector<std::vector<double>> activations_;
        std::vector<std::vector<double>> preactivations_;

        // Construct the time input vector to be appended to the network input.
        // Always includes raw t; if timeEmbeddingDim_ > 0 also appends sinusoidal features
        // [sin(2π·2^i·t), cos(2π·2^i·t)] for i = 0..timeEmbeddingDim_/2-1.
        std::vector<double> timeEmbed(double t) const;

        // Assemble the full network input vector from a (noisy) state vector, the optional
        // conditioning vector, and the diffusion time t. Layout:
        //   [x (dim_), Fourier(x) (Sigma of inputEmbeddingDims_),
        //    cond (conditionDim_), Fourier(cond) (Sigma of conditionEmbeddingDims_),
        //    t (+ time embedding)]
        // Fourier features for coordinate j: [sin(π·2^i·v_j), cos(π·2^i·v_j)] for
        // i = 0..dims[j]/2-1 (none when dims[j] == 0). Depths are per-coordinate.
        // The embedding has no trainable parameters, so back-propagation is unaffected.
        std::vector<double> buildNetworkInput(
            const std::vector<double>& x,
            const std::vector<double>& cond,
            double t
        ) const;

        // Forward pass through the network (training path).
        // Caches activations and pre-activations for the backward pass.
        // input contains state vector, optional condition vector, and diffusion time t.
        std::vector<double> forward(
            const std::vector<double>& x
        );

        // Const forward pass through an arbitrary network (inference path).
        // Does not write to activations_/preactivations_, safe to call during generateSample().
        std::vector<double> forwardInference(
            const std::vector<double>& input,
            const std::vector<Layer>& net
        ) const;

        // Backward pass for gradient computation.
        // Computes gradients w.r.t. network parameters via chain rule.
        //
        // Parameters:
        //   gradOutput - Gradient of loss w.r.t. network output
        //
        // Returns: Gradient of loss w.r.t. network parameters (gradW and gradb for each layer)
        void backward(
            const std::vector<double>& gradOutput
        );

        // Shared reverse-diffusion integrator behind generateSample() and
        // partialReverseSample(): starting from state x at grid time stepStart/steps,
        // iterates the Euler/Heun SDE/ODE updates down to t=0 and returns the
        // de-normalized sample. generateSample() passes stepStart = steps (t=1).
        //
        // Parameters:
        //   x         - Initial normalized state at time stepStart/steps (consumed)
        //   condition - Normalized conditioning vector
        //   stepStart - First reverse step index in [1, steps]; integration covers
        //               t = stepStart/steps, ..., 1/steps
        //   useEMANetworkIfAvailable / useHeun / useSDE / sdeToOdeSigmaThreshold -
        //               identical to generateSample()
        //   steps     - Total number of grid steps (defines dt = 1/steps)
        //
        // Returns: zscore = final state (normalized space), value = de-normalized state
        SBDMGeneratedSample reverseDiffuseFrom(
            std::vector<double> x,
            const std::vector<double>& condition,
            int stepStart,
            bool useEMANetworkIfAvailable,
            bool useHeun,
            bool useSDE,
            int steps,
            double sdeToOdeSigmaThreshold
        );

        // Update network weights using computed gradients (Stochastic Gradient Descent SGD).
        // Applied after backward pass. Alternately, an Adam optimizer step can be implemented
        // in adamUpdate() for better convergence.
        //
        // Parameters:
        //   lr - Learning rate for gradient descent step
        void updateWeights(double lr);

        // Update network weights using Adam optimization algorithm.
        // This method would implement the Adam update rule using the stored gradients and moment estimates.
        //
        // Parameters:
        //   lr - Learning rate for Adam update
        void adamUpdate(double lr);

        // ----- diffusion -----

        // Noise schedule parameter beta(t) over diffusion time [0,1].
        // Linear scheme: interpolation between betaMin_ and betaMax_.
        // Cosine scheme: derived from sigma(t) = sqrt(1 - alpha_bar(t)),
        //     where alpha_bar(t) = cos^2((t + cosineOffset_) / (1 + cosineOffset_) * pi/2).
        //
        // Parameters:
        //   t - Diffusion time parameter in [0,1]
        //
        // Returns: Beta value for the given time step
        double beta(double t) const;

        // Cumulative signal retention factor alpha_bar(t) over diffusion time [0,1].
        // For linear noise schedule, alpha_bar(t) = exp(-integral_0^t beta(s) ds).
        // For cosine noise schedule, alpha_bar(t) = cos^2((t + cosineOffset_) / (1 + cosineOffset_) * pi/2).
        //
        // Parameters:
        //   t - Diffusion time parameter in [0,1]
        //
        // Returns: Cumulative signal retention factor at time t
        double alphabar(double t) const;

        // Cumulative perturbation standard deviation of noise sigma(t) at diffusion time [0,1].
        // Related to the noise schedule via sigma(t) = sqrt(1 - alpha_bar(t)).
        // For LOGSIG: sigma(t) = logSigMin_ * exp(k*t), k = ln(logSigMax_/logSigMin_).
        //
        // Parameters:
        //   t - Diffusion time parameter in [0,1]
        //
        // Returns: Noise standard deviation at time t
        double sigma(double t) const;

        // Time derivative dσ/dt. Used by beta() for the LOGSIG schedule.
        //
        // Parameters:
        //   t - Diffusion time parameter in [0,1]
        //
        // Returns: dσ/dt at time t
        double dSigmadt(double t) const;

        // Add Gaussian noise to a state vector at diffusion time t.
        // Uses external engine (randGaussQ_) for reproducible noise generation.
        // Noise sample is stored in eps for later use in training.
        //
        // Parameters:
        //   x   - Original state vector
        //   t   - Diffusion time parameter in [0,1]
        //   eps - Output: Gaussian noise vector used for perturbation (size = dim_)
        //
        // Returns: Noisy state = sqrt(alpha_bar(t)) * x + sigma(t) * eps
        std::vector<double> addNoise(
            const std::vector<double>& x,
            double t,
            std::vector<double>& eps
        );

        // ----- loss -----

        // Compute Mean Squared Error between predicted score and target score.
        // Used during training to measure model performance.
        //
        // Parameters:
        //   score  - Model's predicted score vector
        //   target - Ground truth target score vector
        //
        // Returns: MSE loss value (scalar)
        double computeLoss(
            const std::vector<double>& score,
            const std::vector<double>& target,
            double weight // use weighted loss to prevent sigma(t) at tiny t from blowing up and dominating the training.
        ) const;

        // Clip gradients to prevent exploding gradients during training.
        //
        // Parameters:
        //   maxNorm - Maximum allowed norm for the gradients. If the total norm exceeds this threshold,
        //             gradients are scaled down to have norm equal to maxNorm.
        void clipGradients(double maxNorm);

        // Apply one EMA step: emaNetwork_ = decay * emaNetwork_ + (1 - decay) * network_.
        // Called after each optimizer step during training.
        void updateEMANetwork();

        // ----- internal vars -----

        // The random engine is NOT owned by this class. It is injected externally
        // by the framework. Below are CLHEP distribution wrappers for actual random number generation.
        // These wrap the engine_ and provide specific probability distributions.
        // - RandFlat:   Uniform distribution on [0,1)
        // - RandGaussQ: Gaussian (normal) distribution with mean=0, sigma=1 (or custom)
        // Both are bound in the constructor to externally managed wrappers.
        CLHEP::RandFlat& randFlat_;       // Used for uniform sampling (e.g., batch selection)
        CLHEP::RandGaussQ& randGaussQ_;   // Used for Gaussian noise in diffusion process

        // Model hyperparameters
        int dim_;               // Dimensionality of state space
        int conditionDim_;      // Dimensionality of the optional conditioning vector
        int timeEmbeddingDim_;  // Sinusoidal time embedding dimensions (0 = raw scalar only)
        std::vector<int> inputEmbeddingDims_;     // Per-coordinate Fourier embedding depth for each state dim (size dim_; 0 = raw only)
        std::vector<int> conditionEmbeddingDims_; // Per-coordinate Fourier embedding depth for each condition dim (size conditionDim_; 0 = raw only)
        int hidden_;            // Hidden layer size
        int layers_;            // Number of network layers

        // Optimizer configuration
        OptimizerType optimizerType_;  // Type of optimizer to use (SGD or ADAM)

        // Adam optimizer parameters
        double adamBeta1_;  // First moment exponential decay rate (default: 0.9)
        double adamBeta2_;  // Second moment exponential decay rate (default: 0.999)
        double adamEps_;    // Small constant for numerical stability (default: 1e-8)

        // Noise schedule configuration
        NoiseScheduleType noiseScheduleType_;  // Type of noise schedule (LINEAR, COSINE, or LOGSIG)

        // Linear noise schedule parameters (beta(t) = betaMin + t*(betaMax - betaMin))
        // default betaMin = 1e-4, betaMax = 0.02 are typical values used in diffusion models, but can be tuned for specific applications.
        double betaMin_;    // Beta value at t=0
        double betaMax_;    // Beta value at t=1

        // Cosine noise schedule parameter (offset to avoid singularity at t=0 in cosine schedule, default: 0.008)
        double cosineOffset_;

        // Log-sigma noise schedule parameters: sigma(t) = logSigMin_ * exp(k*t), k = ln(logSigMax_/logSigMin_)
        double logSigMin_;  // sigma at t=0 (default: 1e-5)
        double logSigMax_;  // sigma at t=1 (default: 1.0)

        // If true, network predicts noise epsilon instead of the score s = -eps/sigma.
        // Prevents value explosion from 1/sigma at small t.
        bool epsPrediction_;

        // Training configuration
        double lossWeightPower_; // Power of the loss function weighting (default: 2.0)
        int batchSize_;  // Batch size for vectorized training (default: 32)
        double gradientClipThreshold_;  // Gradient clipping threshold (default: 1.0)
        double learningRate_; // Learning rate for training (default: 1e-3)

        // Adaptive dimensional weight controller ---
        bool   useDimWeightController_; // if true, per-dimension gradient weights are applied during training
        double dimWeightEMADecay_;       // EMA decay rate for per-dimension loss tracking
        std::vector<double> dimLossEMA_; // per-dim raw MSE EMA (size dim_), updated every epoch
        std::vector<double> dimWeights_; // normalized gradient weights (size dim_), init 1.0

        // EMA copy of network parameters for inference ---
        static constexpr int kEMABatchSizeRef_ = 32; // canonical reference batch size for emaNetworkDecay interpretation
        bool   useEMANetwork_;          // if true, generateSample() uses emaNetwork_ instead of network_
        double emaNetworkDecayBase_;    // user-configured decay at kEMABatchSizeRef_=32 samples/step
        double emaNetworkDecay_;        // effective decay per optimizer step (rescaled from emaNetworkDecayBase_)
        std::vector<Layer> emaNetwork_; // slow-moving EMA copy, only W and b are used

        // In-memory snapshot buffers for the curriculum planner's resume-from-best
        // (see snapshotNetwork()/restoreNetwork()). Not serialized.
        std::vector<Layer> networkSnapshot_;
        std::vector<Layer> emaNetworkSnapshot_;
        std::vector<double> dimWeightsSnapshot_;
        std::vector<double> dimLossEMASnapshot_;
        int   adamStepSnapshot_ = 0;
        bool  hasSnapshot_      = false;

        // Diffusion process discretization
        int diffusionSteps_;  // Number of time steps to generate a sample (default: 200)

        // Training state
        double runningLoss_;  // Accumulated loss for monitoring during training
        int adamStep_; // Step counter for Adam optimizer (used to compute bias-corrected moment estimates)
        size_t trainingSampleSize_; // Total number of training samples

        // Container for tracking training loss over epochs
        std::vector<double> epochLosses_;

        // Variables to track gradient clipping statistics for monitoring
        size_t clipCount_;
        size_t totalClipChecks_;
        double clipScaleAccum_;

        // Data normalization containers
        // Mean and variance of input data (for normalization)
        std::vector<double> dataMean_;
        std::vector<double> dataStdev_;
        // Min and max of normalized data (for possible clamping or rescaling)
        std::vector<double> normMin_;
        std::vector<double> normMax_;
    };
}
