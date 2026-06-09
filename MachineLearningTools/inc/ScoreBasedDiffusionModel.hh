// Module for training and using (variance preserving) score-based diffusion model
// Added by Yongyi Wu
// Mar. 2026

#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
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
            return batchSize_;
        }

        const std::vector<double>& getDimWeights() const { return dimWeights_; }

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
        void train(
            const std::vector<DiffusionTrainingSample>& data,
            int epochs,
            int trainSubsetDataSize = 0,
            bool biasLowSigma = false,
            double tLowBound = 0.0
        );

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
            int diffusionSteps = -1
        );

        // Save the model parameters to a CSV file with annotations for later use.
        // Uses a default filename of "DiffusionModel.csv" if not specified.
        //
        // Parameters:
        //   filename - Path to the CSV file where model parameters will be saved (default: "DiffusionModel.csv")
        void saveModel(const std::string& filename = "DiffusionModel.csv");

        // Load model parameters from a file to restore a previously trained model.
        // If the file contains an [OPTIMIZER_STATE] section (saved by a current saveModel call),
        // the Adam moments and step counter are restored so training can resume seamlessly.
        // Files saved without that section (e.g. older checkpoints) load cleanly with a fresh optimizer state.
        //
        // Parameters:
        //   randFlat / RandGaussQ - CLHEP random number generator wrappers being passed
        //   filename              - Path to the file from which model parameters will be loaded
        static ScoreBasedDiffusionModel loadModel(
            CLHEP::RandFlat& randFlat,
            CLHEP::RandGaussQ& randGaussQ,
            const std::string& filename
        );

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
        bool   useEMANetwork_;      // if true, generateSample() uses emaNetwork_ instead of network_
        double emaNetworkDecay_;    // decay rate applied per optimizer step (e.g. 0.9999)
        std::vector<Layer> emaNetwork_; // slow-moving EMA copy, only W and b are used

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
