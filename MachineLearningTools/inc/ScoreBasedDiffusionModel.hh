// Module for training and using score-based diffusion model
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
            COSINE   // Cosine noise schedule
        };

        // Constructor: Initialize diffusion model with CLHEP random distributions.
        //
        // Parameters:
        //   randFlat              - Reference to CLHEP RandFlat for uniform sampling (externally managed)
        //   randGaussQ            - Reference to CLHEP RandGaussQ for Gaussian noise (externally managed)
        //   dim                   - Dimensionality of the state space
        //   conditionDim          - Dimensionality of the optional conditioning vector (default: 0 for unconditional model)
        //   hidden                - Size of hidden layers in the neural network
        //   layers                - Number of layers in the network
        //   optimizerType         - Type of optimizer to use (SGD or ADAM, default: ADAM)
        //   adamBeta1             - Adam optimizer beta1 parameter (default: 0.9)
        //   adamBeta2             - Adam optimizer beta2 parameter (default: 0.999)
        //   adamEps               - Adam optimizer epsilon parameter (default: 1e-8)
        //   scheduleType          - Type of noise schedule (LINEAR or COSINE, default: COSINE)
        //   betaMin               - Minimum noise schedule parameter (for LINEAR schedule, default: 1e-4)
        //   betaMax               - Maximum noise schedule parameter (for LINEAR schedule, default: 0.02)
        //   cosineOffset          - Offset parameter (for cosine schedule, default: 0.008)
        //   batchSize             - Batch size for training (default: 32)
        //   gradientClipThreshold - Threshold for gradient clipping (default: 1.0)
        //   learningRate          - Learning rate for training (default: 1e-3)
        //   diffusionSteps        - Number of steps in the diffusion process (default: 200)
        ScoreBasedDiffusionModel(
            // Network architecture parameters
            CLHEP::RandFlat& randFlat,
            CLHEP::RandGaussQ& randGaussQ,
            int dim,
            int conditionDim,
            int hidden,
            int layers,
            // Optimizer configuration
            OptimizerType optimizerType = OptimizerType::ADAM,
            double adamBeta1 = 0.9,
            double adamBeta2 = 0.999,
            double adamEps = 1e-8,
            // Noise schedule configuration
            NoiseScheduleType scheduleType = NoiseScheduleType::COSINE,
            double betaMin = 1e-4,
            double betaMax = 0.02,
            double cosineOffset = 0.008,
            // Training configuration
            int batchSize = 32,
            double gradientClipThreshold = 1.0,
            double learningRate = 1e-3,
            // Diffusion process configuration
            int diffusionSteps = 200
        );

        // Train the score network on a batch of samples.
        // Uses random sampling and noise injection via the external engine.
        // Note that training needs to occur on all data samples. Training on multiple small subsets
        // of data and then averaging or aggregating the model parameters is not supported and may lead to 
        // issues as neural networks are not linear. 
        //
        // Parameters:
        //   data         - Training samples (transformed state vectors)
        //   epochs       - Number of training epochs to perform
        void train(
            const std::vector<DiffusionTrainingSample>& data,
            int epochs
        );

        // Generate a new sample from the diffusion model via reverse process.
        // Uses the external random engine for noise generation during sampling.
        //
        // Parameters:
        //   condition       - Optional conditioning vector (must match conditionDim_ when enabled)
        //   useHeun         - If true, uses Heun's method (2nd order). If false, uses Euler method (2nd order, default)
        //   diffusionSteps  - Number of diffusion steps for sampling (default: -1 uses the model's configured diffusionSteps_)
        //
        // Returns: A generated sample vector of dimension dim_
        std::vector<double> generateSample(
            const std::vector<double>& condition = {},
            bool useHeun = true,
            int diffusionSteps = -1
        );

        // Save the model parameters to a CSV file with annotations for later use.
        // Uses a default filename of "DiffusionModel.csv" if not specified.
        //
        // Parameters: 
        //   filename - Path to the CSV file where model parameters will be saved (default: "DiffusionModel.csv")
        void saveModel(const std::string& filename = "DiffusionModel.csv");

        // Load model parameters from a file to restore a previously trained model.
        // Note that as the adam optimizer state is not saved, it is not possible to resume training from a loaded 
        // model and pick up the training process where it left off. The loaded model can only be used for sampling.
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

        // Forward pass through the network.
        // Computes the score function.
        // input actually contains the state vector, optional condition vector, and the time embedding.
        // dimension hence is dim_ + conditionDim_ + 1.
        //
        // Parameters:
        //   x - Input vector of dimension dim_ + conditionDim_ + 1
        //
        // Returns: Output vector (dimension depends on network architecture)
        std::vector<double> forward(
            const std::vector<double>& x
        );

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
        // Linear interpolation between betaMin_ and betaMax_.
        //
        // Parameters:
        //   t - Diffusion time parameter in [0,1]
        //
        // Returns: Beta value for the given time step
        double beta(double t) const;

        // Standard deviation of noise at diffusion time t.
        // Related to the noise schedule via sigma(t) = sqrt(1 - exp(-integral(beta(s) ds))).
        //
        // Parameters:
        //   t - Diffusion time parameter in [0,1]
        //
        // Returns: Noise standard deviation at time t
        double sigma(double t) const;

        // Add Gaussian noise to a state vector at diffusion time t.
        // Uses external engine (randGaussQ_) for reproducible noise generation.
        // Noise sample is stored in eps for later use in training.
        //
        // Parameters:
        //   x   - Original state vector
        //   t   - Diffusion time parameter in [0,1]
        //   eps - Output: Gaussian noise vector used for perturbation (size = dim_)
        //
        // Returns: Noisy state = x + sigma(t) * eps
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
            const std::vector<double>& target
        ) const;

        // Clip gradients to prevent exploding gradients during training.
        //
        // Parameters:
        //   maxNorm - Maximum allowed norm for the gradients. If the total norm exceeds this threshold, 
        //             gradients are scaled down to have norm equal to maxNorm.
        void clipGradients(double maxNorm);

        // ----- internal vars -----
        
        // The random engine is NOT owned by this class. It is injected externally
        // by the framework. Below are CLHEP distribution wrappers for actual random number generation.
        // These wrap the engine_ and provide specific probability distributions.
        // - RandFlat:   Uniform distribution on [0,1)
        // - RandGaussQ: Gaussian (normal) distribution with mean=0, sigma=1 (or custom)
        // Both are initialized in the constructor with the injected engine_.
        CLHEP::RandFlat   randFlat_;      // Used for uniform sampling (e.g., batch selection)
        CLHEP::RandGaussQ randGaussQ_;    // Used for Gaussian noise in diffusion process

        // Model hyperparameters
        int dim_;           // Dimensionality of state space
        int conditionDim_;  // Dimensionality of the optional conditioning vector
        int hidden_;        // Hidden layer size
        int layers_;        // Number of network layers

        // Optimizer configuration
        OptimizerType optimizerType_;  // Type of optimizer to use (SGD or ADAM)

        // Adam optimizer parameters
        double adamBeta1_;  // First moment exponential decay rate (default: 0.9)
        double adamBeta2_;  // Second moment exponential decay rate (default: 0.999)
        double adamEps_;    // Small constant for numerical stability (default: 1e-8)

        // Noise schedule configuration
        NoiseScheduleType noiseScheduleType_;  // Type of noise schedule (LINEAR or COSINE)

        // Linear noise schedule parameters (beta(t) = betaMin + t*(betaMax - betaMin))
        // default betaMin = 1e-4, betaMax = 0.02 are typical values used in diffusion models, but can be tuned for specific applications.
        double betaMin_;    // Beta value at t=0
        double betaMax_;    // Beta value at t=1

        // Cosine noise schedule parameter (offset to avoid singularity at t=0 in cosine schedule, default: 0.008)
        double cosineOffset_;

        // Training configuration
        int batchSize_;  // Batch size for vectorized training (default: 32)
        double gradientClipThreshold_;  // Gradient clipping threshold (default: 1.0)
        double learningRate_; // Learning rate for training (default: 1e-3) 

        // Diffusion process discretization
        int diffusionSteps_;  // Number of time steps to generate a sample (default: 200)

        // Training state
        double runningLoss_;  // Accumulated loss for monitoring during training
        int adamStep_; // Step counter for Adam optimizer (used to compute bias-corrected moment estimates)
        int trainingSampleSize_; // Total number of training samples 

        // Container for tracking training loss over epochs
        std::vector<double> epochLosses_;
        
    };
}
