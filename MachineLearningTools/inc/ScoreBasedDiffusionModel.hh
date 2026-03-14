// Module for training and using score-based diffusion model
// Added by Yongyi Wu
// Mar. 2026

#pragma once

#include <vector>
#include <random>
#include <cmath>
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace mu2e{
    struct DiffusionTrainingSample {
        std::vector<double> x; // transformed state vector (size = dim)
    };

    class ScoreBasedDiffusionModel {
    public:
        // Constructor: Initialize diffusion model with random engine.
        //
        // Parameters:
        //   engine        - Reference to CLHEP random engine (externally managed)
        //   dim           - Dimensionality of the state space
        //   hidden        - Size of hidden layers in the neural network
        //   layers        - Number of layers in the network
        //   betaMin       - Minimum noise schedule parameter
        //   betaMax       - Maximum noise schedule parameter
        //   diffusionSteps - Number of steps in the diffusion process
        ScoreBasedDiffusionModel(
            CLHEP::HepRandomEngine& engine,
            int dim,
            int hidden,
            int layers,
            double betaMin,
            double betaMax,
            int diffusionSteps
        );

        // Train the score network on a batch of samples.
        // Uses random sampling and noise injection via the external engine.
        //
        // Parameters:
        //   data         - Training samples (transformed state vectors)
        //   epochs       - Number of training epochs
        //   batchSize    - Samples per batch (should divide evenly into data.size())
        //   learningRate - Gradient descent step size
        void train(
            const std::vector<DiffusionTrainingSample>& data,
            int epochs,
            int batchSize,
            double learningRate
        );

        // Generate a new sample from the diffusion model via reverse process.
        // Uses the external random engine for noise generation during sampling.
        //
        // Returns: A generated sample vector of dimension dim_
        std::vector<double> generateSample();

    private:

        // ----- network -----

        // Represents a single fully-connected layer with weights and biases.
        struct Layer {
            std::vector<std::vector<double>> W; // Weight matrix [output_size][input_size]
            std::vector<double> b;              // Bias vector [output_size]
        };

        std::vector<Layer> network_; // Network layers in forward order

        // Forward pass through the network.
        // Computes network output given input state.
        //
        // Parameters:
        //   x - Input vector of dimension dim_
        //
        // Returns: Output vector (dimension depends on network architecture)
        std::vector<double> forward(
            const std::vector<double>& x
        );

        // Backward pass for gradient computation.
        // Computes gradients w.r.t. network parameters from output gradients.
        //
        // Parameters:
        //   gradOutput - Gradient of loss w.r.t. network output
        void backward(
            const std::vector<double>& gradOutput
        );

        // Update network weights using computed gradients (SGD).
        // Applied after backward pass.
        //
        // Parameters:
        //   lr - Learning rate for gradient descent step
        void updateWeights(double lr);

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

        // ----- internal vars -----
        // Model hyperparameters
        int dim_;           // Dimensionality of state space
        int hidden_;        // Hidden layer size
        int layers_;        // Number of network layers

        // Noise schedule parameters (beta(t) = betaMin + t*(betaMax - betaMin))
        double betaMin_;    // Beta value at t=0
        double betaMax_;    // Beta value at t=1

        // Diffusion process discretization
        int diffusionSteps_;  // Number of time steps for reverse process

        // The random engine is NOT owned by this class. It is injected externally
        // by the framework. Engine reference must remain valid for the lifetime
        // of this object.
        CLHEP::HepRandomEngine& engine_;

        // CLHEP distribution wrappers for actual random number generation.
        // These wrap the engine_ and provide specific probability distributions.
        // - RandFlat:   Uniform distribution on [0,1)
        // - RandGaussQ: Gaussian (normal) distribution with mean=0, sigma=1 (or custom)
        // Both are initialized in the constructor with the injected engine_.
        CLHEP::RandFlat   randFlat_;      // Used for uniform sampling (e.g., batch selection)
        CLHEP::RandGaussQ randGaussQ_;    // Used for Gaussian noise in diffusion process

        // Training state
        double runningLoss_;  // Accumulated loss for monitoring during training
    };
}
