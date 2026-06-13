#include "Offline/MachineLearningTools/inc/ScoreBasedDiffusionModel.hh"

namespace mu2e {

    // SiLU activation function (Sigmoid Linear Unit)
    // Defined as: silu(x) = x / (1 + exp(-x))
    static double silu(double x) {
        return x / (1.0 + std::exp(-x));
    }

    // Derivative of SiLU activation function, needed for back-propagation.
    // Defined as: silu'(x) = sigmoid(x) * (1 + x * (1 - sigmoid(x)))
    static double siluDeriv(double x) {
        double s = 1.0 / (1.0 + std::exp(-x));
        return s * (1.0 + x * (1.0 - s));
    }

    // Here instantaneous noise scale, cumulative signal retention factor,
    // and cumulative perturbation stddev are defined for variance preserving (VP) SDEs
    // for both linear and cosine noise schedules. This choice is out of the consideration
    // of numerical stability, ease of implementation, and faster sampler convergence
    
    // Instantaneous noise scale (diffusion coefficient) at time t
    double ScoreBasedDiffusionModel::beta(double t) const {
        if (noiseScheduleType_ == NoiseScheduleType::COSINE) {
            // Cosine noise schedule does not use betaMin and betaMax, but we can still define an effective beta if needed
            double f = (t + cosineOffset_) / (1.0 + cosineOffset_);
            double beta = (M_PI / (1.0 + cosineOffset_)) * std::tan(f * M_PI * 0.5); // beta(t) = pi / (1+offset) * tan(pi*f/2)
            // Cap beta to prevent numerical issues at t close to 1
            return std::min(beta, 10.0); // cap beta to a large value
        } else if (noiseScheduleType_ == NoiseScheduleType::LOGSIG) {
            // LOGSIG: beta(t) = 2*sigma*dSigmadt / (1 - sigma^2)
            // Beta diverges as sigma->1 (t->1). For default sigMin=1e-5, k~11.5:
            //   beta=100 is first reached at sigma~0.90 (t~0.97).
            // Capping at 100 keeps beta*dt <= 0.5 for 200 diffusion steps,
            // analogous to the cap=10 used for COSINE (where k is ~7x smaller).
            double s  = sigma(t);
            double sd = dSigmadt(t);
            double denom = 1.0 - s * s;
            if (denom < 1e-12) return 100.0;
            return std::min(2.0 * s * sd / denom, 100.0);
        } else {
            // Linear interpolation of beta(t) between betaMin_ and betaMax_.
            return betaMin_ + t * (betaMax_ - betaMin_);
        }
    }

    // cumulative signal retention factor at time t
    double ScoreBasedDiffusionModel::alphabar(double t) const {
        if (noiseScheduleType_ == NoiseScheduleType::COSINE) {
            // Cosine noise schedule
            double f = (t + cosineOffset_) / (1.0 + cosineOffset_);
            double alpha_bar = std::cos(f * M_PI * 0.5);
            alpha_bar *= alpha_bar; // alpha_bar(t) = cos^2(pi*f/2)
            return alpha_bar;
        } else if (noiseScheduleType_ == NoiseScheduleType::LOGSIG) {
            // LOGSIG: alphabar(t) = 1 - sigma^2(t), clamped to [0,1]
            double s = sigma(t);
            return std::max(0.0, 1.0 - s * s);
        } else {
            // Linear noise schedule
            double integral = betaMin_ * t + 0.5 * (betaMax_ - betaMin_) * t * t;
            return std::exp(-integral); // alpha_bar(t) = exp(-integral(beta(s) ds))
        }
    }

    // Cumulative perturbation stddev (noise level)
    double ScoreBasedDiffusionModel::sigma(double t) const {
        if (noiseScheduleType_ == NoiseScheduleType::LOGSIG) {
            // Direct formula; computed here (not via alphabar) to avoid circular dependency
            double k = std::log(logSigMax_ / logSigMin_);
            return logSigMin_ * std::exp(k * t);
        }
        return std::sqrt(1.0 - alphabar(t)); // sigma(t) = sqrt(1 - alpha_bar(t))
    }

    // Time derivative dσ/dt — only defined for LOGSIG
    double ScoreBasedDiffusionModel::dSigmadt(double t) const {
        if (noiseScheduleType_ != NoiseScheduleType::LOGSIG) {
            throw cet::exception("ScoreBasedDiffusionModel::dSigmadt")
                << "dSigmadt() is only defined for the LOGSIG noise schedule";
        }
        double k = std::log(logSigMax_ / logSigMin_);
        return k * sigma(t);
    }


    std::vector<double> ScoreBasedDiffusionModel::timeEmbed(double t) const {
        std::vector<double> emb;
        emb.reserve(1 + timeEmbeddingDim_);
        emb.push_back(t);
        for (int i = 0; i < timeEmbeddingDim_ / 2; ++i) {
            double freq = 2.0 * M_PI * std::pow(2.0, i);
            emb.push_back(std::sin(freq * t));
            emb.push_back(std::cos(freq * t));
        }
        return emb;
    }

    // Assemble the full network input: raw state coordinates, optional per-coordinate Fourier
    // features, conditioning vector, then the time input (raw t + optional time embedding).
    // The Fourier features give the MLP access to high-frequency structure in the state
    // coordinates (e.g. narrow spectral lines), which a plain MLP is biased against learning.
    std::vector<double> ScoreBasedDiffusionModel::buildNetworkInput(
        const std::vector<double>& x,
        const std::vector<double>& cond,
        double t
    ) const {
        std::vector<double> input;
        input.reserve(dim_ + dim_ * inputEmbeddingDim_
                      + conditionDim_ + conditionDim_ * conditionEmbeddingDim_
                      + 1 + timeEmbeddingDim_);
        // Per-coordinate Fourier features: [sin(π·2^i·v), cos(π·2^i·v)] for i = 0..embDim/2-1
        auto appendFourier = [&input](const std::vector<double>& v, int embDim) {
            for (double vj : v) {
                for (int i = 0; i < embDim / 2; ++i) {
                    double freq = M_PI * std::pow(2.0, i);
                    input.push_back(std::sin(freq * vj));
                    input.push_back(std::cos(freq * vj));
                }
            }
        };
        input.insert(input.end(), x.begin(), x.end());
        if (inputEmbeddingDim_ > 0) appendFourier(x, inputEmbeddingDim_);
        input.insert(input.end(), cond.begin(), cond.end());
        if (conditionEmbeddingDim_ > 0) appendFourier(cond, conditionEmbeddingDim_);
        auto tEmb = timeEmbed(t);
        input.insert(input.end(), tEmb.begin(), tEmb.end());
        return input;
    }

    ScoreBasedDiffusionModel::ScoreBasedDiffusionModel(
        CLHEP::RandFlat& randFlat,
        CLHEP::RandGaussQ& randGaussQ,
        int dim,
        int conditionDim,
        int timeEmbeddingDim,
        int inputEmbeddingDim,
        int conditionEmbeddingDim,
        int hidden,
        int layers,
        OptimizerType optimizerType,
        double adamBeta1,
        double adamBeta2,
        double adamEps,
        NoiseScheduleType scheduleType,
        double betaMin,
        double betaMax,
        double cosineOffset,
        double logSigMin,
        double logSigMax,
        bool epsPrediction,
        double lossWeightPower,
        int batchSize,
        double gradientClipThreshold,
        double learningRate,
        bool useDimWeightController,
        double dimWeightEMADecay,
        bool useEMANetwork,
        double emaNetworkDecay,
        int diffusionSteps,
        bool initializeRandomWeights
    ) : randFlat_(randFlat), randGaussQ_(randGaussQ),
        dim_(dim), conditionDim_(conditionDim), timeEmbeddingDim_(timeEmbeddingDim), inputEmbeddingDim_(inputEmbeddingDim), conditionEmbeddingDim_(conditionEmbeddingDim), hidden_(hidden), layers_(layers),
        optimizerType_(optimizerType),
        adamBeta1_(adamBeta1), adamBeta2_(adamBeta2), adamEps_(adamEps),
        noiseScheduleType_(scheduleType),
        betaMin_(betaMin), betaMax_(betaMax), cosineOffset_(cosineOffset),
        logSigMin_(logSigMin), logSigMax_(logSigMax), epsPrediction_(epsPrediction),
        lossWeightPower_(lossWeightPower), batchSize_(batchSize), gradientClipThreshold_(gradientClipThreshold), learningRate_(learningRate),
        useDimWeightController_(useDimWeightController), dimWeightEMADecay_(dimWeightEMADecay),
        dimLossEMA_(dim, 0.0), dimWeights_(dim, 1.0),
        useEMANetwork_(useEMANetwork), emaNetworkDecayBase_(emaNetworkDecay),
        emaNetworkDecay_(useEMANetwork ? std::pow(emaNetworkDecay, (double)batchSize / kEMABatchSizeRef_) : emaNetworkDecay),
        diffusionSteps_(diffusionSteps),
        runningLoss_(0.0), adamStep_(0), trainingSampleSize_(0), epochLosses_(),
        clipCount_(0), totalClipChecks_(0), clipScaleAccum_(0.0),
        dataMean_(dim + conditionDim, 0.0), dataStdev_(dim + conditionDim, 1.0),
        normMin_(dim + conditionDim, -999.0), normMax_(dim + conditionDim, 999.0) {

        // Validate model dimensions and parameters
        if (dim <= 0 || conditionDim < 0 || hidden <= 0 || layers <= 0) {
            throw cet::exception("ScoreBasedDiffusionModel::initialization") << "Invalid model dimensions";
        }
        if (timeEmbeddingDim_ != 0 && (timeEmbeddingDim_ < 2 || timeEmbeddingDim_ % 2 != 0)) {
            throw cet::exception("ScoreBasedDiffusionModel::initialization")
                << "timeEmbeddingDim must be 0 (raw scalar) or an even integer >= 2";
        }
        if (inputEmbeddingDim_ != 0 && (inputEmbeddingDim_ < 2 || inputEmbeddingDim_ % 2 != 0)) {
            throw cet::exception("ScoreBasedDiffusionModel::initialization")
                << "inputEmbeddingDim must be 0 (raw coordinates) or an even integer >= 2";
        }
        if (conditionEmbeddingDim_ != 0 && (conditionEmbeddingDim_ < 2 || conditionEmbeddingDim_ % 2 != 0)) {
            throw cet::exception("ScoreBasedDiffusionModel::initialization")
                << "conditionEmbeddingDim must be 0 (raw values) or an even integer >= 2";
        }
        if (conditionEmbeddingDim_ > 0 && conditionDim_ == 0) {
            throw cet::exception("ScoreBasedDiffusionModel::initialization")
                << "conditionEmbeddingDim > 0 requires conditionDim > 0";
        }
        if (batchSize <= 0) {
            throw cet::exception("ScoreBasedDiffusionModel::initialization") << "Invalid batchSize";
        }
        if (diffusionSteps <= 0) {
            throw cet::exception("ScoreBasedDiffusionModel::initialization") << "Invalid diffusionSteps";
        }

        // ------------------------------------------------------------
        // Network architecture
        //
        // Input dimension = dim_ (+ per-coordinate Fourier features) + conditionDim_ + 1
        // (+1 because diffusion time t is appended to the input vector)
        // ------------------------------------------------------------

        int inputSize = dim_ + dim_ * inputEmbeddingDim_
                      + conditionDim_ + conditionDim_ * conditionEmbeddingDim_
                      + 1 + timeEmbeddingDim_;
        int in = inputSize;

        // Weight initialization scale (local constant so it can be tuned easily)
        // const double weightInitScale = 0.02; // not scalable with size, can lead to instability for larger models and too slow training for smaller models
        // const double weightInitScale = std::sqrt(2.0 / in); // He initialization for ReLU activations, SiLU is similar to ReLU in terms of variance preservation
        double reducedHe = 0.5 * std::sqrt(2.0 / in); // Scaled He initialization found to improve stability
        const double weightInitScale = std::min(reducedHe, 0.3); // Cap the weight initialization scale to prevent instability for very small input sizes

        for (int l = 0; l < layers_; ++l) {

            // Last layer outputs the score vector (dimension = dim_)
            int out = (l == layers_ - 1) ? dim_ : hidden_;

            Layer layer;

            // Allocate weights
            layer.W.resize(out, std::vector<double>(in, 0.0));

            // Allocate biases
            layer.b.resize(out, 0.0);

            // Allocate gradient buffers
            layer.gradW.resize(out, std::vector<double>(in, 0.0));
            layer.gradb.resize(out, 0.0);

            // Allocate Adam buffers
            layer.mW.resize(out, std::vector<double>(in, 0.0));
            layer.vW.resize(out, std::vector<double>(in, 0.0));

            layer.mb.resize(out, 0.0);
            layer.vb.resize(out, 0.0);

            if (initializeRandomWeights) {
                // --------------------------------------------------------
                // Weight initialization
                //
                // Small Gaussian initialization improves training stability
                // --------------------------------------------------------
                for (int i = 0; i < out; ++i) {
                    for (int j = 0; j < in; ++j) {
                        layer.W[i][j] = weightInitScale * randGaussQ_.fire();
                        // this generates a Gaussian random number with mean=0 and sigma=weightInitScale
                    }
                    layer.b[i] = 0.0;
                }
            }

            if (layer.W.empty() || layer.W[0].size() != static_cast<size_t>(in) || //Check that weight matrix has correct input dimension
                layer.W.size() != static_cast<size_t>(out) ||                      //Check that weight matrix has correct output dimension
                layer.b.size() != static_cast<size_t>(out)) {                      //Check that bias vector has correct dimension
                throw cet::exception("ScoreBasedDiffusionModel::initialization")
                    << "Layer shape initialization mismatch";
            }

            network_.push_back(layer);

            // Output becomes next layer input
            in = out;
        }

        // Initialize emaNetwork_ with the same W/b as network_ so it starts from the same point
        emaNetwork_.resize(network_.size());
        for (size_t l = 0; l < network_.size(); ++l) {
            emaNetwork_[l].W = network_[l].W;
            emaNetwork_[l].b = network_[l].b;
        }

        // Print layer and model configuration
        std::ostringstream oss;
        oss << "ScoreBasedDiffusionModel initialized\n"
            << "Model configuration:\n"
            << "  Network architecture:\n"
            << "    | dim=" << dim_ << "\n"
            << "    | conditionDim=" << conditionDim_ << "\n"
            << "    | timeEmbeddingDim=" << timeEmbeddingDim_ << (timeEmbeddingDim_ == 0 ? " (raw scalar)" : " (sinusoidal)") << "\n"
            << "    | inputEmbeddingDim=" << inputEmbeddingDim_ << (inputEmbeddingDim_ == 0 ? " (raw coordinates)" : " (Fourier features per state coordinate)") << "\n"
            << "    | conditionEmbeddingDim=" << conditionEmbeddingDim_ << (conditionEmbeddingDim_ == 0 ? " (raw values)" : " (Fourier features per condition coordinate)") << "\n"
            << "    | hidden=" << hidden_ << "\n"
            << "    | layers=" << layers_ << "\n"
            << "  Optimizer configuration:\n"
            << "    | Optimizer=" << (optimizerType_ == OptimizerType::ADAM ? "Adam" : "SGD") << "\n";
        if (optimizerType_ == OptimizerType::ADAM) {
            oss << "      |- AdamBeta1=" << adamBeta1_ << "\n"
                << "      |- AdamBeta2=" << adamBeta2_ << "\n"
                << "      |- AdamEps=" << adamEps_ << "\n";
        }
        oss << "  Noise schedule configuration:\n";
        if (noiseScheduleType_ == NoiseScheduleType::COSINE) {
            oss << "    | NoiseSchedule=Cosine\n"
                << "      |- CosineOffset=" << cosineOffset_ << "\n";
        } else if (noiseScheduleType_ == NoiseScheduleType::LOGSIG) {
            oss << "    | NoiseSchedule=LogSig\n"
                << "      |- LogSigMin=" << logSigMin_ << "\n"
                << "      |- LogSigMax=" << logSigMax_ << "\n";
        } else {
            oss << "    | NoiseSchedule=Linear\n"
                << "      |- BetaMin=" << betaMin_ << "\n"
                << "      |- BetaMax=" << betaMax_ << "\n";
        }
        oss << "  Training configuration:\n"
            << "    | EpsPrediction=" << (epsPrediction_ ? "true" : "false") << "\n"
            << "    | LossWeightPower=" << lossWeightPower_ << "\n"
            << "    | BatchSize=" << batchSize_ << "\n"
            << "    | GradientClipThreshold=" << gradientClipThreshold_ << "\n"
            << "    | LearningRate=" << learningRate_ << "\n"
            << "    | DimWeightController=" << (useDimWeightController_ ? "enabled" : "disabled");
        if (useDimWeightController_)
            oss << " (EMADecay=" << dimWeightEMADecay_ << ")";
        oss << "\n"
            << "  EMA network:\n"
            << "    | EMANetwork=" << (useEMANetwork_ ? "enabled" : "disabled");
        if (useEMANetwork_)
            oss << " (decay=" << emaNetworkDecay_ << ")";
        oss << "\n"
            << "  Diffusion process configuration:\n"
            << "    | DiffusionSteps=" << diffusionSteps_ << "\n";
        mf::LogInfo("ScoreBasedDiffusionModel::initialize") << oss.str();

        // Reserve space for forward-pass caches
        activations_.reserve(layers_ + 1);
        preactivations_.reserve(layers_);
    }

    // normalize input data
    // note the order in the stdev and mean is always x then cond
    void ScoreBasedDiffusionModel::normalizeData(
        const std::vector<double>& mean,
        const std::vector<double>& stdev,
        std::vector<DiffusionTrainingSample>& data
    )
    {
        if (data.empty()) {
            throw cet::exception("ScoreBasedDiffusionModel::normalizeData")
                << "No training data provided";
        }

        const size_t totalDim = dim_ + conditionDim_;

        // check mean and stdev dimensions
        if (mean.size() != totalDim || stdev.size() != totalDim) {
            throw cet::exception("ScoreBasedDiffusionModel::normalizeData")
                << "Mean or standard deviation dimension mismatch, expected "
                << totalDim << " but got " << mean.size()
                << " for mean and " << stdev.size() << " for standard deviation";
        }
        // store the mean and standard deviation of the original data
        dataMean_ = mean;
        dataStdev_ = stdev;

        // containers to track the mean and standard deviation of the normalized data
        std::vector<double> normMean(totalDim, 0.0);
        std::vector<double> normStdev(totalDim, 0.0);
        // containers to track the sum of squares of the normalized data
        std::vector<double> M2(totalDim, 0.0);
        // reset min/max trackers
        normMin_.assign(totalDim, std::numeric_limits<double>::max());
        normMax_.assign(totalDim, std::numeric_limits<double>::lowest());

        size_t count = 0;

        // Normalize the training data using the given mean and standard deviation
        // keep track of the mean and standard deviation, min and max values
        // of the normalized data
        for (auto& sample : data) {
            ++count;

            // deal with sample.x first
            for (int i = 0; i < dim_; ++i) {
                if (stdev[i] == 0.0) {
                    throw cet::exception("ScoreBasedDiffusionModel::normalizeData")
                        << "Zero standard deviation encountered at dimension "
                        << i;
                }

                double value = (sample.x[i] - mean[i]) / stdev[i];
                sample.x[i] = value; //overwrite the original value

                // update min/max values
                normMin_[i] = std::min(normMin_[i], value);
                normMax_[i] = std::max(normMax_[i], value);
                // Welford update
                double delta = value - normMean[i];
                normMean[i] += delta / static_cast<double>(count);
                double delta2 = value - normMean[i];
                M2[i] += delta * delta2;
            }
            //now deal with sample.cond
            for (int i = 0; i < conditionDim_; ++i) {
                const int idx = dim_ + i;
                if (stdev[idx] == 0.0) {
                    throw cet::exception("ScoreBasedDiffusionModel::normalizeData")
                        << "Zero standard deviation encountered at dimension "
                        << idx;
                }

                double value = (sample.cond[i] - mean[idx]) / stdev[idx];
                sample.cond[i] = value;

                // update min/max
                normMin_[idx] = std::min(normMin_[idx], value);
                normMax_[idx] = std::max(normMax_[idx], value);
                // Welford update
                double delta = value - normMean[idx];
                normMean[idx] += delta / static_cast<double>(count);
                double delta2 = value - normMean[idx];
                M2[idx] += delta * delta2;
            }
        }

        for (size_t i = 0; i < totalDim; ++i) {
            normStdev[i] = std::sqrt(M2[i] / static_cast<double>(count));
            mf::LogInfo("ScoreBasedDiffusionModel::normalizeData")
                << "Dimension " << i << " normalized: mean = " << normMean[i] << ", stdev = " << normStdev[i];

            // sanity checks
            if (std::abs(normMean[i]) > 0.1) {
                mf::LogWarning("ScoreBasedDiffusionModel::normalizeData")
                    << "Normalized mean deviates from 0 at dimension "
                    << i
                    << ": mean = " << normMean[i];
            }
            if (std::abs(normStdev[i] - 1.0) > 0.1) {
                mf::LogWarning("ScoreBasedDiffusionModel::normalizeData")
                    << "Normalized stdev deviates from 1 at dimension "
                    << i
                    << ": stdev = " << normStdev[i];
            }
        }
        return;
    }

    // forward pass to compute the score function given input (state + time embedding)
    std::vector<double> ScoreBasedDiffusionModel::forward(
        const std::vector<double>& input
    )
    {
        if (network_.empty() || network_[0].W.empty() || network_[0].W[0].empty()) {
            throw cet::exception("ScoreBasedDiffusionModel::forward")
                << "Network is not properly initialized";
        }
        if (input.size() != network_[0].W[0].size()) {
            throw cet::exception("ScoreBasedDiffusionModel::forward")
                << "Input dimension mismatch: got " << input.size()
                << ", expected " << network_[0].W[0].size();
        }

        activations_.clear();
        preactivations_.clear();

        // Input vector contains both the state and the time embedding (e.g., concatenated together).
        std::vector<double> x = input;

        // Cache the input as the activation of layer 0 for use in back-propagation.
        activations_.push_back(x);

        // Forward pass through the network layers
        for (size_t l = 0; l < network_.size(); ++l)
        {
            // Compute pre-activation z = W*x + b for the current layer
            auto& layer = network_[l];
            std::vector<double> z(layer.W.size());
            for (size_t i = 0; i < layer.W.size(); ++i)
            {
                double v = layer.b[i];
                for (size_t j = 0; j < layer.W[i].size(); ++j)
                {
                    v += layer.W[i][j]*x[j];
                }
                z[i] = v;
            }

            // Cache pre-activation for back-propagation
            preactivations_.push_back(z);

            // Apply activation function (SiLU) to get the output of the current layer.
            std::vector<double> y(z.size());
            if (l != network_.size() - 1) { // No activation on the last layer (output layer), as it predicts the score directly.
                for (size_t i = 0; i < z.size(); ++i) {
                    y[i] = silu(z[i]);
                }
            }
            else
                y = z;

            // Cache activation for back-propagation
            activations_.push_back(y);

            // Output of current layer becomes input to the next layer
            x = y;
        }

        return x;
    }

    // Const forward pass for inference — does not cache activations or preactivations.
    // Used by generateSample() so it can use either network_ or emaNetwork_ without
    // corrupting the activations needed by the training backward pass.
    std::vector<double> ScoreBasedDiffusionModel::forwardInference(
        const std::vector<double>& input,
        const std::vector<Layer>& net
    ) const
    {
        if (net.empty() || net[0].W.empty() || net[0].W[0].empty()) {
            throw cet::exception("ScoreBasedDiffusionModel::forwardInference")
                << "Network is not properly initialized";
        }
        if (input.size() != net[0].W[0].size()) {
            throw cet::exception("ScoreBasedDiffusionModel::forwardInference")
                << "Input dimension mismatch: got " << input.size()
                << ", expected " << net[0].W[0].size();
        }

        std::vector<double> x = input;
        for (size_t l = 0; l < net.size(); ++l) {
            const auto& layer = net[l];
            std::vector<double> z(layer.W.size());
            for (size_t i = 0; i < layer.W.size(); ++i) {
                double v = layer.b[i];
                for (size_t j = 0; j < layer.W[i].size(); ++j)
                    v += layer.W[i][j] * x[j];
                z[i] = v;
            }
            if (l != net.size() - 1) {
                for (size_t i = 0; i < z.size(); ++i)
                    z[i] = silu(z[i]);
            }
            x = z;
        }
        return x;
    }

    // backward pass to compute gradients of the loss w.r.t. network parameters using chain rule
    void ScoreBasedDiffusionModel::backward(
        const std::vector<double>& gradOutput
    )
    {
        std::vector<double> grad = gradOutput;

        // Back-propagate through layers in reverse order
        for (int l = static_cast<int>(network_.size()) - 1; l >= 0; l--)
        {
            auto& layer = network_[l];

            auto& aPrev = activations_[l];
            auto& z = preactivations_[l];

            std::vector<double> gradZ(grad.size());

            // For the output layer, the gradient is directly from the loss w.r.t. output (no activation function).
            if (l == static_cast<int>(network_.size()) - 1) {
                gradZ = grad;
            } else { // For hidden layers, apply the derivative of the activation function (SiLU) to the gradient.
                for (size_t i = 0; i < grad.size(); ++i) {
                    gradZ[i] = grad[i] * siluDeriv(z[i]);
                }
            }

            // Compute gradients w.r.t. weights and biases for the current layer
            for (size_t i = 0; i < layer.W.size(); i++) {
                for (size_t j = 0; j < layer.W[i].size(); ++j) {
                    // Gradient of loss w.r.t. weight W[i][j]
                    layer.gradW[i][j] += gradZ[i] * aPrev[j];
                }

                // Gradient of loss w.r.t. bias b[i]
                // Note: We accumulate gradients in gradW and gradb because we will apply the optimizer step
                // after processing a batch of samples.
                layer.gradb[i] += gradZ[i];
            }

            if (layer.W.empty() || layer.W[0].empty()) {
                throw cet::exception("ScoreBasedDiffusionModel::backward")
                    << "Encountered empty layer weights during backward pass";
            }
            // Compute gradient w.r.t. input of the current layer to propagate to the previous layer
            std::vector<double> gradPrev(layer.W[0].size(), 0.0);

            for (size_t j = 0; j < gradPrev.size(); ++j) {
                for (size_t i = 0; i < layer.W.size(); ++i) {
                    gradPrev[j] += layer.W[i][j] * gradZ[i];
                }
            }

            grad = gradPrev;
        }
    }

    // Update network weights using computed gradients (Stochastic Gradient Descent).
    void ScoreBasedDiffusionModel::updateWeights(
        double lr
    )
    {
        for (auto& layer : network_)
        {
            if (layer.W.empty() || layer.W[0].empty() || layer.b.empty()) {
                throw cet::exception("ScoreBasedDiffusionModel::updateWeights")
                    << "Encountered malformed layer during SGD update";
            }

            size_t outSize = layer.W.size();
            size_t inSize  = layer.W[0].size();

            // Update weights
            for (size_t i = 0; i < outSize; ++i)
            {
                for (size_t j = 0; j < inSize; ++j)
                {
                    // Update weight using SGD: W = W - lr * gradW
                    layer.W[i][j] -= lr * layer.gradW[i][j];

                    // reset gradient
                    layer.gradW[i][j] = 0.0;
                }

                // Update bias
                layer.b[i] -= lr * layer.gradb[i];

                // reset gradient
                layer.gradb[i] = 0.0;
            }
        }
    }

    // Update network weights using Adam optimization algorithm.
    // Adam optimizer maintains running estimates of the first and second moments of the gradients (m and v) for each parameter,
    // and uses these to adapt the learning rate for each parameter individually.
    // The update rule for a parameter theta with gradient g at time step t is:
    // m_t = b1 * m_{t-1} + (1 - b1) * g
    // v_t = b2 * v_{t-1} + (1 - b2) * g^2
    // mhat = m_t / (1 - b1^t)  (bias-corrected first moment estimate)
    // vhat = v_t / (1 - b2^t)  (bias-corrected second moment estimate)
    // theta = theta - lr * mhat / (sqrt(vhat) + eps)
    void ScoreBasedDiffusionModel::adamUpdate(
        double lr
    )
    {
        adamStep_++; // initialize to 0 in constructor,
                     // increment first to avoid division by zero in the first step when computing bias-corrected estimates
                     // (1 - b1^t) and (1 - b2^t). Do not modify this to increment after the update.

        for (auto& layer : network_)
        {
            if (layer.W.empty() || layer.W[0].empty() || layer.b.empty()) {
                throw cet::exception("ScoreBasedDiffusionModel::adamUpdate")
                    << "Encountered malformed layer during Adam update";
            }

            for (size_t i = 0; i < layer.W.size(); ++i) {
                // Update weights using Adam update rule
                for (size_t j = 0; j < layer.W[i].size(); ++j) {
                    double g = layer.gradW[i][j];

                    layer.mW[i][j] = adamBeta1_ * layer.mW[i][j] + (1 - adamBeta1_) * g;
                    layer.vW[i][j] = adamBeta2_ * layer.vW[i][j] + (1 - adamBeta2_) * g * g;

                    double mhat = layer.mW[i][j] / (1 - std::pow(adamBeta1_, adamStep_));
                    double vhat = layer.vW[i][j] / (1 - std::pow(adamBeta2_, adamStep_));

                    layer.W[i][j] -= lr * mhat / (std::sqrt(vhat) + adamEps_);

                    layer.gradW[i][j] = 0.0;
                }

                // Update bias using Adam update rule
                double g = layer.gradb[i];

                layer.mb[i] = adamBeta1_ * layer.mb[i] + (1 - adamBeta1_) * g;
                layer.vb[i] = adamBeta2_ * layer.vb[i] + (1 - adamBeta2_) * g * g;

                double mhat = layer.mb[i] / (1 - std::pow(adamBeta1_, adamStep_));
                double vhat = layer.vb[i] / (1 - std::pow(adamBeta2_, adamStep_));

                layer.b[i] -= lr * mhat / (std::sqrt(vhat) + adamEps_);

                layer.gradb[i] = 0.0;
            }
        }
    }

    // Add noise to the input state x at diffusion time t to produce the noisy sample xt.
    std::vector<double> ScoreBasedDiffusionModel::addNoise(
        const std::vector<double>& x,
        double t, // Diffusion time parameter in [0,1]
        std::vector<double>& eps
    ){
        if (x.size() != static_cast<size_t>(dim_)) {
            throw cet::exception("ScoreBasedDiffusionModel::addNoise")
                << "Input x dimension mismatch: got " << x.size()
                << ", expected " << dim_;
        }

        double alpha_bar, s;
        if (noiseScheduleType_ == NoiseScheduleType::LOGSIG) {
            s= sigma(t);
            alpha_bar = 1.0 - s * s; // directly compute alpha_bar from sigma for efficiency
        } else {
            alpha_bar = alphabar(t);
            s = std::sqrt(1.0 - alpha_bar); // i.e. sigma(t), avoid recalculating alpha_bar inside the function for efficiency
        }

        // Generate Gaussian noise vector eps of dimension dim_ using the external random engine.
        eps.resize(dim_);

        std::vector<double> xt(dim_);

        for (int i = 0; i < dim_; ++i) {
            eps[i] = randGaussQ_.fire(); // Gaussian N(0,1)
            xt[i] = std::sqrt(alpha_bar) * x[i] + s * eps[i];
        }

        return xt;
    }

    // Compute Mean Squared Error per dimension between predicted score and target score.
    double ScoreBasedDiffusionModel::computeLoss(
        const std::vector<double>& score,
        const std::vector<double>& target,
        double weight // use weighted loss to prevent sigma(t) at tiny t from blowing up and dominating the training.
    ) const {

        if (score.size() != static_cast<size_t>(dim_)) {
            throw cet::exception("ScoreBasedDiffusionModel::computeLoss")
                << "Score dimension mismatch: got " << score.size()
                << ", expected " << dim_;
        }
        if (target.size() != static_cast<size_t>(dim_)) {
            throw cet::exception("ScoreBasedDiffusionModel::computeLoss")
                << "Target dimension mismatch: got " << target.size()
                << ", expected " << dim_;
        }

        double loss = 0.0;

        for (int i = 0; i < dim_; ++i) {
            double d = score[i] - target[i];
            loss += weight * d * d;
        }

        return loss / dim_;
    }

    // Clip gradients to prevent exploding gradients during training. This is done by scaling down
    // the gradients if their L2 norm exceeds a specified threshold.
    void ScoreBasedDiffusionModel::clipGradients(
        double maxNorm
    )
    {
        double norm = 0.0;

        for (auto& layer : network_)
        {
            for (auto& row : layer.gradW)
                for (double g : row)
                    norm += g*g;
            for (double g : layer.gradb)
                norm += g*g;
        }
        norm = std::sqrt(norm);

        totalClipChecks_++;
        // If the norm is below the threshold, no clipping is needed.
        if (norm <= maxNorm) {
            clipScaleAccum_ += 1.0;
            return;
        }
        clipCount_++;

        // Scale down the gradients to have the specified maximum norm.
        double scale = maxNorm / norm;
        clipScaleAccum_ += scale;
        for (auto& layer : network_)
        {
            for (auto& row : layer.gradW)
                for (auto& g : row)
                    g *= scale;
            for (auto& g : layer.gradb)
                g *= scale;
        }
    }

    // Train the diffusion model using the provided training data. The training loop iterates over epochs and batches,
    // applying noise to the input samples, computing the score predictions, calculating the loss, and performing
    // back-propagation to update the model parameters.
    // For a data sample of 5M 6-dimensional vectors of double precision (8 Byte), total memory for the data is 5M*6*8 = 240 MB,
    // which is manageable for in-memory training.
    // Much larger datasets may require streaming from disk or using mini-batches that do not fit entirely in memory.
    void ScoreBasedDiffusionModel::updateEMANetwork()
    {
        for (size_t l = 0; l < network_.size(); ++l) {
            for (size_t i = 0; i < network_[l].W.size(); ++i) {
                for (size_t j = 0; j < network_[l].W[i].size(); ++j)
                    emaNetwork_[l].W[i][j] = emaNetworkDecay_ * emaNetwork_[l].W[i][j]
                                           + (1.0 - emaNetworkDecay_) * network_[l].W[i][j];
                emaNetwork_[l].b[i] = emaNetworkDecay_ * emaNetwork_[l].b[i]
                                    + (1.0 - emaNetworkDecay_) * network_[l].b[i];
            }
        }
    }

    void ScoreBasedDiffusionModel::train(
        const std::vector<DiffusionTrainingSample>& data,
        int epochs,
        int trainSubsetDataSize,
        bool biasLowSigma,
        double tLowBound,
        double tFocusLow,
        double tFocusHigh,
        double tFocusFraction
    )
    {
        // Check that the network is properly initialized before training.
        if (network_.empty() || network_[0].W.empty() || network_[0].W[0].empty()) {
            throw cet::exception("ScoreBasedDiffusionModel::train")
                << "Network is not properly initialized";
        }

        if (tLowBound >= 1.0 || tLowBound < 0.0) {
            throw cet::exception("ScoreBasedDiffusionModel::train")
                << "Invalid tLowBound: must be in [0,1)";
        }

        if (tFocusFraction < 0.0 || tFocusFraction > 1.0) {
            throw cet::exception("ScoreBasedDiffusionModel::train")
                << "Invalid tFocusFraction: must be in [0,1]";
        }
        if (tFocusFraction > 0.0 &&
            (tFocusLow < 0.0 || tFocusHigh > 1.0 || tFocusLow >= tFocusHigh)) {
            throw cet::exception("ScoreBasedDiffusionModel::train")
                << "Invalid t focus window: require 0 <= tFocusLow < tFocusHigh <= 1 (got ["
                << tFocusLow << ", " << tFocusHigh << "])";
        }

        const double eps_safe = 1e-12;
        const size_t N = data.size();

        trainingSampleSize_ = N; // Store the total number of training samples

        if (N <= 1)
        {
            throw cet::exception("ScoreBasedDiffusionModel::train")
                << "Insufficient data samples for training (samples must be > 1)";
        }

        for (size_t i = 0; i < N; ++i) {
            if (data[i].x.size() != static_cast<size_t>(dim_)) {
                throw cet::exception("ScoreBasedDiffusionModel::train")
                    << "Training sample " << i << " has x dimension " << data[i].x.size()
                    << ", expected " << dim_;
            }
            if (data[i].cond.size() != static_cast<size_t>(conditionDim_)) {
                throw cet::exception("ScoreBasedDiffusionModel::train")
                    << "Training sample " << i << " has conditioning dimension " << data[i].cond.size()
                    << ", expected " << conditionDim_;
            }
        }

        // Determine how many samples to use per epoch after shuffling
        const size_t activeN = (trainSubsetDataSize > 0 && static_cast<size_t>(trainSubsetDataSize) < N)
                               ? static_cast<size_t>(trainSubsetDataSize)
                               : N;

        // Create index vector once
        std::vector<size_t> indices(N);
        for (size_t i = 0; i < N; ++i)
            indices[i] = i;

        // Data shuffling and batching is performed at the epoch level to ensure that each epoch sees the data
        // in a different order, which can improve training convergence.
        for (int e = 0; e < epochs; ++e)
        {
            // Shuffle indices using Fisher–Yates and RandFlat for reproducibility
            for (size_t i = N - 1; i > 0; --i)
            {
                size_t j = static_cast<size_t>(randFlat_.fire() * (i + 1));
                j = std::min(j, i); // Ensure j is within bounds
                std::swap(indices[i], indices[j]);
            }

            // Counter for number of samples processed in the epoch (used for averaging loss). Avoid using N directly to
            // allow for early stopping or partial epoch processing if needed.
            int n = 0;

            int batchCounter = 0;
            double epochLoss = 0.0;
            std::vector<double> epochDimLoss(dim_, 0.0);

            // iterate over the shuffled data samples (subset if trainSubsetDataSize was specified)
            for (size_t idx = 0; idx < activeN; ++idx)
            {
                const auto& sample = data[indices[idx]];

                // Sample diffusion time. If tFocusFraction > 0, this sample draws t uniformly from the
                // focus window [tFocusLow, tFocusHigh] with probability tFocusFraction, concentrating
                // gradient steps on a target sigma band while the remaining samples keep full coverage.
                // Otherwise: if tLowBound is set, we scale the uniform random number to be in
                // [tLowBound, 1] instead of [0,1] to focus training on later diffusion steps with higher
                // noise levels, which can improve stability and convergence in some cases. If biasLowSigma
                // is true, we apply a quadratic transformation to the sampled time to bias it towards smaller
                // values (less noise), which can help the model learn the score function more accurately at
                // early diffusion steps where the signal is stronger. Not meant to be used together.

                double t = randFlat_.fire();
                if (tFocusFraction > 0.0 && randFlat_.fire() < tFocusFraction) {
                    t = tFocusLow + (tFocusHigh - tFocusLow) * t; // scale t to the focus window
                } else {
                    if (tLowBound > 0.0) {
                        t = tLowBound + (1.0 - tLowBound) * t; // scale t to [tLowBound, 1]
                    }
                    if (biasLowSigma) {
                        t = t * t; // focus on smaller t value (smaller sigma)
                    }
                }
                // container for noise vector eps
                std::vector<double> eps;

                auto xt = addNoise(sample.x,t,eps);
                double s = std::max(sigma(t), eps_safe); // add eps_safe to prevent division by zero in case of very small sigma
                double weight = std::pow(s, lossWeightPower_); // sigma(t)^power as the weight. Quadratic power good for convergence but missing details

                // Target: either the noise epsilon (epsPrediction_=true) or the score -eps/sigma(t).
                // Epsilon prediction prevents 1/sigma blow-up at small t.
                std::vector<double> target(dim_);
                for (int i = 0; i < dim_; ++i) {
                    target[i] = epsPrediction_ ? eps[i] : -eps[i] / s;
                }

                // Prepare input vector for the network from the noisy sample xt, condition, and time.
                auto input = buildNetworkInput(xt, sample.cond, t);

                // Check that input dimension matches the network's expected input dimension
                if (input.size() != network_[0].W[0].size()) {
                    throw cet::exception("ScoreBasedDiffusionModel::train")
                        << "Training input dimension mismatch: got " << input.size()
                        << ", expected " << network_[0].W[0].size();
                }

                // Forward pass to compute the predicted score from the network given the input (noisy sample + time).
                auto score = forward(input);

                // Compute the loss (Mean Squared Error) per dimension between the predicted score and the target score.
                double loss = computeLoss(score, target, weight);

                // Compute the gradient of the loss w.r.t. the predicted score.
                // dimWeights_[i] rescales each dimension's gradient; normalized to mean=1 so overall
                // gradient scale is preserved. Accumulate raw squared residuals for the controller EMA.
                std::vector<double> grad(dim_);
                for (int i = 0; i < dim_; ++i) {
                    double residual = score[i] - target[i];
                    if (useDimWeightController_)
                        epochDimLoss[i] += residual * residual;
                    grad[i] = 2.0 * weight * dimWeights_[i] * residual / dim_;
                }

                // Backward pass to compute gradients of the loss w.r.t. network parameters using the computed gradient of the loss w.r.t. the predicted score.
                backward(grad);

                // Increment batch counter and apply optimizer step if batch size is reached. This allows for vectorized updates after processing a batch of
                // samples, which can improve training efficiency and convergence.
                batchCounter++;
                if (batchCounter == batchSize_) {
                    // Average accumulated gradients so optimizer step size is independent of batch size.
                    const double invBatch = 1.0 / static_cast<double>(batchCounter);
                    for (auto& layer : network_) {
                        for (auto& row : layer.gradW)
                            for (auto& g : row)
                                g *= invBatch;
                        for (auto& g : layer.gradb)
                            g *= invBatch;
                    }

                    clipGradients(gradientClipThreshold_); // Clip gradients to prevent exploding gradients

                    // Apply optimizer based on configuration
                    if (optimizerType_ == OptimizerType::ADAM) {
                        adamUpdate(learningRate_);
                    } else {
                        updateWeights(learningRate_);
                    }
                    if (useEMANetwork_) updateEMANetwork();
                    batchCounter = 0;
                }

                epochLoss += loss; // Accumulate loss for monitoring.
                n++;
            }

            // Final flush: apply optimizer to remaining gradients
            if (batchCounter > 0) {
                // Average by the true final mini-batch size (which may be < batchSize_).
                const double invBatch = 1.0 / static_cast<double>(batchCounter);
                for (auto& layer : network_) {
                    for (auto& row : layer.gradW)
                        for (auto& g : row)
                            g *= invBatch;
                    for (auto& g : layer.gradb)
                        g *= invBatch;
                }

                clipGradients(gradientClipThreshold_); // Clip gradients before final update

                if (optimizerType_ == OptimizerType::ADAM) {
                    adamUpdate(learningRate_);
                } else {
                    updateWeights(learningRate_);
                }
                if (useEMANetwork_) updateEMANetwork();
            }

            // Update per-dimension EMA loss and recompute gradient weights
            if (useDimWeightController_ && n > 0) {
                double emaSum = 0.0;
                for (int i = 0; i < dim_; ++i) {
                    dimLossEMA_[i] = dimWeightEMADecay_ * dimLossEMA_[i]
                                   + (1.0 - dimWeightEMADecay_) * (epochDimLoss[i] / n);
                    emaSum += dimLossEMA_[i];
                }
                double emaMean = emaSum / dim_;
                if (emaMean > 0.0) {
                    for (int i = 0; i < dim_; ++i)
                        dimWeights_[i] = dimLossEMA_[i] / emaMean;
                }
                std::ostringstream woss;
                woss << "Epoch " << e << " dimWeights: [";
                for (int i = 0; i < dim_; ++i) {
                    woss << dimWeights_[i];
                    if (i < dim_ - 1) woss << ", ";
                }
                woss << "]";
                mf::LogInfo("ScoreBasedDiffusionModel::train") << woss.str();
            }

            epochLoss /= n;
            double clipRatio = (totalClipChecks_ > 0) ? double(clipCount_) / totalClipChecks_ : 0.0;
            double avgClipScale = (totalClipChecks_ > 0) ? clipScaleAccum_ / totalClipChecks_ : 1.0;

            mf::LogInfo("ScoreBasedDiffusionModel::train") << "Epoch " << e << "  Loss=" << epochLoss << "  ClipRatio=" << clipRatio << "  AvgClipScale=" << avgClipScale;
            epochLosses_.push_back(epochLoss);
        }
    }

    void ScoreBasedDiffusionModel::saveModel(
        const std::string& filename
    )
    {
        // Write binary
        std::ofstream out(filename, std::ios::binary);
        if (!out) {
            throw cet::exception("ScoreBasedDiffusionModel::saveModel")
                << "Cannot open file " << filename;
        }

        auto wI32 = [&](int32_t  v){ out.write(reinterpret_cast<const char*>(&v), 4); };
        auto wU32 = [&](uint32_t v){ out.write(reinterpret_cast<const char*>(&v), 4); };
        auto wU64 = [&](uint64_t v){ out.write(reinterpret_cast<const char*>(&v), 8); };
        auto wF64 = [&](double   v){ out.write(reinterpret_cast<const char*>(&v), 8); };
        auto wVec = [&](const std::vector<double>& v){
            wU64(static_cast<uint64_t>(v.size()));
            if (!v.empty()) out.write(reinterpret_cast<const char*>(v.data()),
                                      static_cast<std::streamsize>(v.size() * sizeof(double)));
        };
        auto wMat = [&](const std::vector<std::vector<double>>& m){
            // Writes a matrix as outSize rows, each row preceded by its column count.
            // Row dimensions may differ (though in practice they don't).
            wU64(static_cast<uint64_t>(m.size()));
            for (const auto& row : m) wVec(row);
        };

        // Magic + version
        // Version history: 1 = original format; 2 = adds inputEmbeddingDim and
        // conditionEmbeddingDim after timeEmbeddingDim.
        out.write("SBDM", 4);
        wU32(2u);

        // Model architecture & hyper-parameters
        wI32(static_cast<int32_t>(dim_));
        wI32(static_cast<int32_t>(conditionDim_));
        wI32(static_cast<int32_t>(timeEmbeddingDim_));
        wI32(static_cast<int32_t>(inputEmbeddingDim_));
        wI32(static_cast<int32_t>(conditionEmbeddingDim_));
        wI32(static_cast<int32_t>(hidden_));
        wI32(static_cast<int32_t>(layers_));
        wI32(optimizerType_ == OptimizerType::ADAM ? 0 : 1);
        wF64(adamBeta1_); wF64(adamBeta2_); wF64(adamEps_);
        int32_t schedIdx = (noiseScheduleType_ == NoiseScheduleType::LINEAR) ? 0
                         : (noiseScheduleType_ == NoiseScheduleType::COSINE)  ? 1 : 2;
        wI32(schedIdx);
        wF64(betaMin_); wF64(betaMax_); wF64(cosineOffset_);
        wF64(logSigMin_); wF64(logSigMax_);
        wI32(epsPrediction_ ? 1 : 0);
        wF64(lossWeightPower_);
        wI32(static_cast<int32_t>(batchSize_));
        wF64(gradientClipThreshold_);
        wF64(learningRate_);
        wI32(useDimWeightController_ ? 1 : 0);
        wF64(dimWeightEMADecay_);
        wI32(useEMANetwork_ ? 1 : 0);
        wF64(emaNetworkDecay_);
        wI32(static_cast<int32_t>(diffusionSteps_));

        // Network weights
        wU32(static_cast<uint32_t>(network_.size()));
        for (const auto& layer : network_) {
            wU32(static_cast<uint32_t>(layer.W.size()));
            wU32(static_cast<uint32_t>(layer.W[0].size()));
            for (const auto& row : layer.W)
                out.write(reinterpret_cast<const char*>(row.data()),
                          static_cast<std::streamsize>(row.size() * sizeof(double)));
            out.write(reinterpret_cast<const char*>(layer.b.data()),
                      static_cast<std::streamsize>(layer.b.size() * sizeof(double)));
        }

        // Data normalisation
        wVec(dataMean_); wVec(dataStdev_); wVec(normMin_); wVec(normMax_);

        // Training history
        wU64(static_cast<uint64_t>(epochLosses_.size()));
        wU64(static_cast<uint64_t>(trainingSampleSize_));
        for (double v : epochLosses_) wF64(v);

        // Optimizer state
        wI32(static_cast<int32_t>(adamStep_));
        for (const auto& layer : network_) {
            wMat(layer.mW); wMat(layer.vW);
            wVec(layer.mb); wVec(layer.vb);
        }
        wVec(dimLossEMA_); wVec(dimWeights_);

        // EMA network (optional section flagged by a boolean)
        wI32(useEMANetwork_ ? 1 : 0);
        if (useEMANetwork_) {
            wU32(static_cast<uint32_t>(emaNetwork_.size()));
            for (const auto& layer : emaNetwork_) {
                wU32(static_cast<uint32_t>(layer.W.size()));
                wU32(static_cast<uint32_t>(layer.W[0].size()));
                for (const auto& row : layer.W)
                    out.write(reinterpret_cast<const char*>(row.data()),
                              static_cast<std::streamsize>(row.size() * sizeof(double)));
                out.write(reinterpret_cast<const char*>(layer.b.data()),
                          static_cast<std::streamsize>(layer.b.size() * sizeof(double)));
            }
        }
    }

    void ScoreBasedDiffusionModel::saveModelCsv(
        const std::string& filename
    )
    {
        std::ofstream out(filename);

        if (!out) {
            throw cet::exception("ScoreBasedDiffusionModel::saveModelCsv")
                << "Cannot open file " << filename;
        }

        // Write CSV header and model configuration parameters with annotations
        out << "[MODEL_PARAMETERS]\n";
        out << "Parameter,Value\n";
        // Model architecture
        out << "dim," << dim_ << "\n";
        out << "conditionDim," << conditionDim_ << "\n";
        out << "timeEmbeddingDim," << timeEmbeddingDim_ << "\n";
        out << "inputEmbeddingDim," << inputEmbeddingDim_ << "\n";
        out << "conditionEmbeddingDim," << conditionEmbeddingDim_ << "\n";
        out << "hidden," << hidden_ << "\n";
        out << "layers," << layers_ << "\n";
        // Optimizer configuration
        out << "optimizerType," << (optimizerType_ == OptimizerType::ADAM ? "ADAM" : "SGD") << "\n";
        out << "adamBeta1," << adamBeta1_ << "\n";
        out << "adamBeta2," << adamBeta2_ << "\n";
        out << "adamEps," << adamEps_ << "\n";
        // Noise schedule configuration
        std::string schedName;
        if (noiseScheduleType_ == NoiseScheduleType::COSINE) schedName = "COSINE";
        else if (noiseScheduleType_ == NoiseScheduleType::LOGSIG) schedName = "LOGSIG";
        else schedName = "LINEAR";
        out << "noiseScheduleType," << schedName << "\n";
        out << "betaMin," << betaMin_ << "\n";
        out << "betaMax," << betaMax_ << "\n";
        out << "cosineOffset," << cosineOffset_ << "\n";
        out << "logSigMin," << logSigMin_ << "\n";
        out << "logSigMax," << logSigMax_ << "\n";
        // Training configuration
        out << "epsPrediction," << (epsPrediction_ ? "1" : "0") << "\n";
        out << "lossWeightPower," << lossWeightPower_ << "\n";
        out << "batchSize," << batchSize_ << "\n";
        out << "gradientClipThreshold," << gradientClipThreshold_ << "\n";
        out << "learningRate," << learningRate_ << "\n";
        // Dimensional weight controller
        out << "useDimWeightController," << (useDimWeightController_ ? "1" : "0") << "\n";
        out << "dimWeightEMADecay," << dimWeightEMADecay_ << "\n";
        // EMA network
        out << "useEMANetwork," << (useEMANetwork_ ? "1" : "0") << "\n";
        out << "emaNetworkDecay," << emaNetworkDecay_ << "\n";
        // Diffusion process configuration
        out << "diffusionSteps," << diffusionSteps_ << "\n";

        // Write network architecture header
        out << "\n[NETWORK_PARAMETERS]\n";
        out << "numLayers," << network_.size() << "\n";
        out << std::fixed << std::setprecision(17); // Use fixed-point notation with high precision for weights and biases

        for (size_t layerIdx = 0; layerIdx < network_.size(); ++layerIdx) {
            // Write layer dimensions
            auto& layer = network_[layerIdx];
            size_t outSize = layer.W.size();
            size_t inSize  = layer.W[0].size();
            out << "\nLayer" << layerIdx << "_OutSize," << outSize << "\n";
            out << "Layer" << layerIdx << "_InSize," << inSize << "\n";

            // Write weights in matrix format
            out << "Layer" << layerIdx << "_Weights\n";
            for (const auto& row : layer.W) {
                for (size_t j = 0; j < row.size(); ++j) {
                    out << row[j];
                    if (j < row.size() - 1) out << ",";
                }
                out << "\n";
            }

            // Write biases
            out << "Layer" << layerIdx << "_Biases\n";
            for (size_t j = 0; j < layer.b.size(); ++j) {
                out << layer.b[j];
                if (j < layer.b.size() - 1) out << ",";
            }
            out << "\n";
        }

        // Write data normalization parameters
        out << "\n[DATA_NORMALIZATION]\n";
        out << "numDimensions," << dataMean_.size() << "\n";
        out << "dataMean\n";
        for (size_t i = 0; i < dataMean_.size(); ++i) {
            out << dataMean_[i];
            if (i < dataMean_.size() - 1) out << ",";
        }
        out << "\n";
        out << "dataStdev\n";
        for (size_t i = 0; i < dataStdev_.size(); ++i) {
            out << dataStdev_[i];
            if (i < dataStdev_.size() - 1) out << ",";
        }
        out << "\n";
        out << "normMin\n";
        for (size_t i = 0; i < normMin_.size(); ++i) {
            out << normMin_[i];
            if (i < normMin_.size() - 1) out << ",";
        }
        out << "\n";
        out << "normMax\n";
        for (size_t i = 0; i < normMax_.size(); ++i) {
            out << normMax_[i];
            if (i < normMax_.size() - 1) out << ",";
        }
        out << "\n";

        // Write training history
        out << "\n[TRAINING_HISTORY]\n";
        out << "numEpochs," << epochLosses_.size() << "\n";
        out << "trainingSampleSize," << trainingSampleSize_ << "\n";
        out << "EpochNumber,Loss\n";
        for (size_t i = 0; i < epochLosses_.size(); ++i) {
            out << i << "," << epochLosses_[i] << "\n";
        }

        // Write Adam optimizer state so training can be resumed from this checkpoint
        out << "\n[OPTIMIZER_STATE]\n";
        out << "adamStep," << adamStep_ << "\n";
        for (size_t layerIdx = 0; layerIdx < network_.size(); ++layerIdx) {
            auto& layer = network_[layerIdx];

            out << "\nLayer" << layerIdx << "_mW\n";
            for (const auto& row : layer.mW) {
                for (size_t j = 0; j < row.size(); ++j) {
                    out << row[j];
                    if (j < row.size() - 1) out << ",";
                }
                out << "\n";
            }

            out << "Layer" << layerIdx << "_vW\n";
            for (const auto& row : layer.vW) {
                for (size_t j = 0; j < row.size(); ++j) {
                    out << row[j];
                    if (j < row.size() - 1) out << ",";
                }
                out << "\n";
            }

            out << "Layer" << layerIdx << "_mb\n";
            for (size_t j = 0; j < layer.mb.size(); ++j) {
                out << layer.mb[j];
                if (j < layer.mb.size() - 1) out << ",";
            }
            out << "\n";

            out << "Layer" << layerIdx << "_vb\n";
            for (size_t j = 0; j < layer.vb.size(); ++j) {
                out << layer.vb[j];
                if (j < layer.vb.size() - 1) out << ",";
            }
            out << "\n";
        }

        // Dim weight controller state
        out << "dimLossEMA\n";
        for (int i = 0; i < dim_; ++i) {
            out << dimLossEMA_[i];
            if (i < dim_ - 1) out << ",";
        }
        out << "\n";
        out << "dimWeights\n";
        for (int i = 0; i < dim_; ++i) {
            out << dimWeights_[i];
            if (i < dim_ - 1) out << ",";
        }
        out << "\n";

        // EMA network weights (used for inference)
        if (useEMANetwork_) {
            out << "\n[EMA_NETWORK]\n";
            out << "numLayers," << emaNetwork_.size() << "\n";
            for (size_t layerIdx = 0; layerIdx < emaNetwork_.size(); ++layerIdx) {
                auto& layer = emaNetwork_[layerIdx];
                out << "\nLayer" << layerIdx << "_OutSize," << layer.W.size() << "\n";
                out << "Layer" << layerIdx << "_InSize," << layer.W[0].size() << "\n";

                out << "Layer" << layerIdx << "_Weights\n";
                for (const auto& row : layer.W) {
                    for (size_t j = 0; j < row.size(); ++j) {
                        out << row[j];
                        if (j < row.size() - 1) out << ",";
                    }
                    out << "\n";
                }

                out << "Layer" << layerIdx << "_Biases\n";
                for (size_t j = 0; j < layer.b.size(); ++j) {
                    out << layer.b[j];
                    if (j < layer.b.size() - 1) out << ",";
                }
                out << "\n";
            }
        }
    }

    ScoreBasedDiffusionModel ScoreBasedDiffusionModel::loadModel(
        CLHEP::RandFlat& randFlat,
        CLHEP::RandGaussQ& randGaussQ,
        const std::string& filename
    )
    {
        // Dispatch based on file extension
        bool isBin = filename.size() >= 4 &&
                     filename.compare(filename.size() - 4, 4, ".bin") == 0;
        bool isCsv = filename.size() >= 4 &&
                     filename.compare(filename.size() - 4, 4, ".csv") == 0;
        if (isBin) {
            std::ifstream bin(filename, std::ios::binary);
            if (!bin) {
                throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                    << "Cannot open binary file " << filename;
            }

            // Fail with an attributable message on a short read, instead of letting
            // an uninitialized/garbage value propagate into a vector allocation (which
            // surfaces as an opaque std::length_error / std::bad_alloc).
            auto fail = [&](const char* ctx) {
                throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                    << "Truncated or corrupt binary checkpoint " << filename
                    << ": unexpected EOF or read error while reading " << ctx << ".";
            };
            auto rI32 = [&]() -> int32_t  { int32_t  v = 0; bin.read(reinterpret_cast<char*>(&v), 4); if (!bin) fail("int32");  return v; };
            auto rU32 = [&]() -> uint32_t { uint32_t v = 0; bin.read(reinterpret_cast<char*>(&v), 4); if (!bin) fail("uint32"); return v; };
            auto rU64 = [&]() -> uint64_t { uint64_t v = 0; bin.read(reinterpret_cast<char*>(&v), 8); if (!bin) fail("uint64"); return v; };
            auto rF64 = [&]() -> double   { double   v = 0; bin.read(reinterpret_cast<char*>(&v), 8); if (!bin) fail("double"); return v; };
            // Sanity-bound a size/count field read from the file BEFORE it is used to size
            // a container. A misaligned-but-not-yet-EOF stream yields a self-consistent
            // header followed by garbage counts; bounding them here converts that into a
            // clear error rather than an allocation that throws length_error.
            constexpr uint64_t kMaxBinElems = 100000000ull; // 1e8 doubles (~800 MB) per dimension
            auto checkCount = [&](uint64_t v, const char* ctx) -> uint64_t {
                if (v > kMaxBinElems)
                    throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                        << "Corrupt binary checkpoint " << filename << ": implausible " << ctx
                        << " count " << v << " (max " << kMaxBinElems
                        << "). File is likely truncated or damaged.";
                return v;
            };
            auto rVec = [&](size_t n) -> std::vector<double> {
                uint64_t stored = rU64();
                if (stored != n)
                    throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                        << "Binary vector size mismatch: expected " << n << ", got " << stored;
                std::vector<double> v(n);
                if (n) {
                    bin.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(double)));
                    if (!bin) fail("vector payload");
                }
                return v;
            };
            auto rMat = [&](size_t rows, size_t cols) -> std::vector<std::vector<double>> {
                uint64_t storedRows = rU64();
                if (storedRows != rows)
                    throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                        << "Binary matrix row count mismatch: expected " << rows << ", got " << storedRows;
                std::vector<std::vector<double>> m(rows);
                for (auto& row : m) {
                    uint64_t storedCols = rU64();
                    if (storedCols != cols)
                        throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                            << "Binary matrix col count mismatch: expected " << cols << ", got " << storedCols;
                    row.resize(cols);
                    bin.read(reinterpret_cast<char*>(row.data()), static_cast<std::streamsize>(cols * sizeof(double)));
                    if (!bin) fail("matrix payload");
                }
                return m;
            };

            // Magic + version
            char magic[4]; bin.read(magic, 4);
            if (std::string(magic, 4) != "SBDM")
                throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                    << "Invalid magic bytes in binary file " << filename;
            uint32_t version = rU32();
            if (version != 1 && version != 2)
                throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                    << "Unsupported binary file version " << version;

            // Model parameters
            int dim            = rI32();
            int conditionDim   = rI32();
            int timeEmbeddingDim = rI32();
            // Version 2 adds the Fourier input embeddings; version-1 files predate them
            int inputEmbeddingDim     = (version >= 2) ? rI32() : 0;
            int conditionEmbeddingDim = (version >= 2) ? rI32() : 0;
            int hidden         = rI32();
            int layers         = rI32();
            int optType        = rI32();
            OptimizerType optimizerType = (optType == 0) ? OptimizerType::ADAM : OptimizerType::SGD;
            double adamBeta1   = rF64(), adamBeta2 = rF64(), adamEps = rF64();
            int schedIdx       = rI32();
            NoiseScheduleType scheduleType = (schedIdx == 0) ? NoiseScheduleType::LINEAR
                                           : (schedIdx == 1) ? NoiseScheduleType::COSINE
                                                             : NoiseScheduleType::LOGSIG;
            double betaMin     = rF64(), betaMax = rF64(), cosineOffset = rF64();
            double logSigMin   = rF64(), logSigMax = rF64();
            bool epsPrediction = (rI32() != 0);
            double lossWeightPower = rF64();
            int batchSize      = rI32();
            double gradientClipThreshold = rF64();
            double learningRate = rF64();
            bool useDimWeightController = (rI32() != 0);
            double dimWeightEMADecay = rF64();
            bool useEMANetwork = (rI32() != 0);
            double emaNetworkDecay = rF64();
            int diffusionSteps = rI32();

            // Network weights
            uint32_t numLayers = static_cast<uint32_t>(checkCount(rU32(), "network layer"));
            std::vector<Layer> loadedNetwork(numLayers);
            for (auto& layer : loadedNetwork) {
                uint32_t outSize = static_cast<uint32_t>(checkCount(rU32(), "layer outSize"));
                uint32_t inSize  = static_cast<uint32_t>(checkCount(rU32(), "layer inSize"));
                layer.W.resize(outSize, std::vector<double>(inSize));
                for (auto& row : layer.W) {
                    bin.read(reinterpret_cast<char*>(row.data()),
                             static_cast<std::streamsize>(inSize * sizeof(double)));
                    if (!bin) fail("network weight row");
                }
                layer.b.resize(outSize);
                bin.read(reinterpret_cast<char*>(layer.b.data()),
                         static_cast<std::streamsize>(outSize * sizeof(double)));
                if (!bin) fail("network biases");
            }

            // Data normalisation
            uint64_t normDim = checkCount(rU64(), "normalisation dimension");
            std::vector<double> dataMean(normDim), dataStdev(normDim), normMin(normDim), normMax(normDim);
            bin.read(reinterpret_cast<char*>(dataMean.data()),  static_cast<std::streamsize>(normDim * sizeof(double)));
            bin.read(reinterpret_cast<char*>(dataStdev.data()), static_cast<std::streamsize>(normDim * sizeof(double)));
            bin.read(reinterpret_cast<char*>(normMin.data()),   static_cast<std::streamsize>(normDim * sizeof(double)));
            bin.read(reinterpret_cast<char*>(normMax.data()),   static_cast<std::streamsize>(normDim * sizeof(double)));
            if (!bin) fail("normalisation arrays");

            // Training history
            uint64_t numEpochs        = checkCount(rU64(), "epoch-loss");
            uint64_t trainingSampleSz = rU64();
            std::vector<double> epochLosses(numEpochs);
            for (auto& v : epochLosses) v = rF64();

            // Optimizer state
            int loadedAdamStep = rI32();
            std::vector<std::vector<std::vector<double>>> loadedMW(numLayers), loadedVW(numLayers);
            std::vector<std::vector<double>> loadedMb(numLayers), loadedVb(numLayers);
            for (size_t l = 0; l < numLayers; ++l) {
                size_t outSize = loadedNetwork[l].W.size();
                size_t inSize  = loadedNetwork[l].W[0].size();
                loadedMW[l] = rMat(outSize, inSize);
                loadedVW[l] = rMat(outSize, inSize);
                loadedMb[l] = rVec(outSize);
                loadedVb[l] = rVec(outSize);
            }
            std::vector<double> loadedDimLossEMA = rVec(static_cast<size_t>(dim));
            std::vector<double> loadedDimWeights  = rVec(static_cast<size_t>(dim));

            // EMA network
            bool hasEMA = (rI32() != 0);
            std::vector<Layer> loadedEmaNetwork;
            if (hasEMA) {
                uint32_t emaLayers = static_cast<uint32_t>(checkCount(rU32(), "EMA layer"));
                loadedEmaNetwork.resize(emaLayers);
                for (auto& layer : loadedEmaNetwork) {
                    uint32_t outSize = static_cast<uint32_t>(checkCount(rU32(), "EMA layer outSize"));
                    uint32_t inSize  = static_cast<uint32_t>(checkCount(rU32(), "EMA layer inSize"));
                    layer.W.resize(outSize, std::vector<double>(inSize));
                    for (auto& row : layer.W) {
                        bin.read(reinterpret_cast<char*>(row.data()),
                                 static_cast<std::streamsize>(inSize * sizeof(double)));
                        if (!bin) fail("EMA weight row");
                    }
                    layer.b.resize(outSize);
                    bin.read(reinterpret_cast<char*>(layer.b.data()),
                             static_cast<std::streamsize>(outSize * sizeof(double)));
                    if (!bin) fail("EMA biases");
                }
            }

            if (!bin) {
                throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                    << "Read error or unexpected EOF in binary file " << filename;
            }

            // Validate
            if (dim <= 0 || conditionDim < 0 || hidden <= 0 || layers <= 0)
                throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                    << "Invalid model parameters in binary file";
            if ((int)normDim != dim + conditionDim)
                throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                    << "Normalisation dimension mismatch in binary file";

            // Reconstruct model
            ScoreBasedDiffusionModel model(
                randFlat, randGaussQ, dim, conditionDim, timeEmbeddingDim,
                inputEmbeddingDim, conditionEmbeddingDim, hidden, layers,
                optimizerType, adamBeta1, adamBeta2, adamEps, scheduleType,
                betaMin, betaMax, cosineOffset, logSigMin, logSigMax,
                epsPrediction, lossWeightPower, batchSize, gradientClipThreshold, learningRate,
                useDimWeightController, dimWeightEMADecay, useEMANetwork, emaNetworkDecay,
                diffusionSteps, false
            );

            for (size_t l = 0; l < model.network_.size(); ++l) {
                model.network_[l].W = loadedNetwork[l].W;
                model.network_[l].b = loadedNetwork[l].b;
            }
            model.dataMean_  = dataMean;
            model.dataStdev_ = dataStdev;
            model.normMin_   = normMin;
            model.normMax_   = normMax;
            model.epochLosses_       = epochLosses;
            model.trainingSampleSize_ = static_cast<size_t>(trainingSampleSz);

            // Restore optimizer state
            model.adamStep_ = loadedAdamStep;
            for (size_t l = 0; l < model.network_.size(); ++l) {
                model.network_[l].mW = loadedMW[l];
                model.network_[l].vW = loadedVW[l];
                model.network_[l].mb = loadedMb[l];
                model.network_[l].vb = loadedVb[l];
            }
            model.dimLossEMA_ = loadedDimLossEMA;
            model.dimWeights_ = loadedDimWeights;

            // Restore EMA network
            if (hasEMA && useEMANetwork) {
                for (size_t l = 0; l < model.emaNetwork_.size(); ++l) {
                    model.emaNetwork_[l].W = loadedEmaNetwork[l].W;
                    model.emaNetwork_[l].b = loadedEmaNetwork[l].b;
                }
                mf::LogInfo("ScoreBasedDiffusionModel::loadModel")
                    << "EMA network weights restored from binary checkpoint.";
            } else if (useEMANetwork) {
                for (size_t l = 0; l < model.emaNetwork_.size(); ++l) {
                    model.emaNetwork_[l].W = model.network_[l].W;
                    model.emaNetwork_[l].b = model.network_[l].b;
                }
                mf::LogInfo("ScoreBasedDiffusionModel::loadModel")
                    << "No EMA network in binary checkpoint — initialized from network weights.";
            }

            mf::LogInfo("ScoreBasedDiffusionModel::loadModel")
                << "Binary model loaded from " << filename
                << " (adamStep=" << loadedAdamStep << ", epochs=" << numEpochs << ")";
            return model;
        } else if (isCsv) {

            std::ifstream in(filename);

            if (!in) {
                throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                    << "Cannot open file " << filename;
            }

            std::string line;

            // Temporary storage
            int dim = 0, conditionDim = 0, timeEmbeddingDim = 0, hidden = 0, layers = 0;
            int inputEmbeddingDim = 0, conditionEmbeddingDim = 0; // absent in files saved before the Fourier input embedding
            OptimizerType optimizerType = OptimizerType::ADAM;
            double adamBeta1 = 0.0, adamBeta2 = 0.0, adamEps = 0.0;
            NoiseScheduleType scheduleType = NoiseScheduleType::COSINE;
            double betaMin = 0.0, betaMax = 0.0, cosineOffset = 0.0;
            double logSigMin = 1e-5, logSigMax = 1.0;
            bool epsPrediction = false;
            double lossWeightPower = 2.0;
            int batchSize = 1, diffusionSteps = 1;
            double gradientClipThreshold = 0.0, learningRate = 0.0;

            std::vector<double> dataMean, dataStdev, normMin, normMax;

            std::vector<Layer> loadedNetwork;
            std::vector<double> epochLosses;

            // Optimizer state (optional — absent in files saved before this feature was added)
            int loadedAdamStep = 0;
            bool optimizerStateLoaded = false;
            std::vector<std::vector<std::vector<double>>> loadedMW, loadedVW;
            std::vector<std::vector<double>> loadedMb, loadedVb;
            std::vector<double> loadedDimLossEMA, loadedDimWeights;

            // New controller / EMA network parameters (default to constructor defaults if absent)
            bool useDimWeightController = false;
            double dimWeightEMADecay = 0.99;
            bool useEMANetwork = true;
            double emaNetworkDecay = 0.9999;

            // EMA network weights (optional)
            std::vector<Layer> loadedEmaNetwork;
            bool emaNetworkLoaded = false;

            // Helper lambda to split CSV
            auto split = [](const std::string& s) {
                std::vector<std::string> tokens;
                std::stringstream ss(s);
                std::string item;
                while (std::getline(ss, item, ',')) {
                    // trim whitespace from item (beginning and end, should not be present if csv is generated by code, but just in case)
                    item.erase(item.begin(), std::find_if(item.begin(), item.end(), [](unsigned char ch) { return !std::isspace(ch); }));
                    item.erase(std::find_if(item.rbegin(), item.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), item.end());
                    tokens.push_back(item);
                }
                return tokens;
            };
            // Helper lambda to extract layer index from strings like "Layer10_OutSize"
            auto getLayerIdx = [](const std::string& s) {
                size_t start = 5; // after "Layer"
                size_t end = s.find('_', start);
                if (end == std::string::npos)
                    throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Malformed layer string: " << s;
                return std::stoi(s.substr(start, end - start));
            };

            // Parse file
            std::string section;

            while (std::getline(in, line)) {
                if (line.empty()) continue;

                // Detect section headers
                if (line[0] == '[') {
                    section = line;
                    continue;
                }

                // MODEL PARAMETERS
                if (section == "[MODEL_PARAMETERS]") {
                    if (line == "Parameter,Value") continue;

                    auto tokens = split(line);
                    if (tokens.size() != 2) throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "ScoreBasedDiffusionModel::loadModel: invalid line in [MODEL_PARAMETERS] section: " << line;

                    const std::string& key = tokens[0];
                    const std::string& val = tokens[1];

                    if (key == "dim") dim = std::stoi(val);
                    else if (key == "conditionDim") conditionDim = std::stoi(val);
                    else if (key == "timeEmbeddingDim") timeEmbeddingDim = std::stoi(val);
                    else if (key == "inputEmbeddingDim") inputEmbeddingDim = std::stoi(val);
                    else if (key == "conditionEmbeddingDim") conditionEmbeddingDim = std::stoi(val);
                    else if (key == "hidden") hidden = std::stoi(val);
                    else if (key == "layers") layers = std::stoi(val);
                    else if (key == "optimizerType")
                        optimizerType = (val == "ADAM") ? OptimizerType::ADAM : OptimizerType::SGD;
                    else if (key == "adamBeta1") adamBeta1 = std::stod(val);
                    else if (key == "adamBeta2") adamBeta2 = std::stod(val);
                    else if (key == "adamEps") adamEps = std::stod(val);
                    else if (key == "noiseScheduleType") {
                        if (val == "COSINE") scheduleType = NoiseScheduleType::COSINE;
                        else if (val == "LOGSIG") scheduleType = NoiseScheduleType::LOGSIG;
                        else scheduleType = NoiseScheduleType::LINEAR;
                    }
                    else if (key == "betaMin") betaMin = std::stod(val);
                    else if (key == "betaMax") betaMax = std::stod(val);
                    else if (key == "cosineOffset") cosineOffset = std::stod(val);
                    else if (key == "logSigMin") logSigMin = std::stod(val);
                    else if (key == "logSigMax") logSigMax = std::stod(val);
                    else if (key == "epsPrediction") epsPrediction = (val == "1");
                    else if (key == "lossWeightPower") lossWeightPower = std::stod(val);
                    else if (key == "batchSize") batchSize = std::stoi(val);
                    else if (key == "gradientClipThreshold") gradientClipThreshold = std::stod(val);
                    else if (key == "learningRate") learningRate = std::stod(val);
                    else if (key == "useDimWeightController") useDimWeightController = (val == "1");
                    else if (key == "dimWeightEMADecay") dimWeightEMADecay = std::stod(val);
                    else if (key == "useEMANetwork") useEMANetwork = (val == "1");
                    else if (key == "emaNetworkDecay") emaNetworkDecay = std::stod(val);
                    else if (key == "diffusionSteps") diffusionSteps = std::stoi(val);
                }

                // NETWORK PARAMETERS
                else if (section == "[NETWORK_PARAMETERS]") {

                    if (line.rfind("numLayers", 0) == 0) {
                        int numLayers = std::stoi(split(line)[1]);
                        loadedNetwork.resize(numLayers);
                        if (loadedNetwork.size() != (size_t)layers) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Layer count mismatch";
                        }
                        continue;
                    }

                    // Layer sizes
                    if (line.find("_OutSize") != std::string::npos) {
                        auto tokens = split(line);
                        int outSize = std::stoi(tokens[1]);
                        int layerIdx = getLayerIdx(tokens[0]);
                        loadedNetwork[layerIdx].W.resize(outSize);
                        loadedNetwork[layerIdx].b.resize(outSize);
                    }
                    else if (line.find("_InSize") != std::string::npos) {
                        auto tokens = split(line);
                        int inSize = std::stoi(tokens[1]);
                        int layerIdx = getLayerIdx(tokens[0]);
                        int expectedInSize = dim + dim * inputEmbeddingDim
                                           + conditionDim + conditionDim * conditionEmbeddingDim
                                           + 1 + timeEmbeddingDim;
                        if (layerIdx == 0 && inSize != expectedInSize) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                                << "Input layer input size mismatch (expected "
                                << expectedInSize << ", got " << inSize << ")";
                        }
                        if (loadedNetwork[layerIdx].W.empty()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "InSize encountered before OutSize for layer " << layerIdx;
                        }
                        for (auto& row : loadedNetwork[layerIdx].W)
                            row.resize(inSize);
                    }
                    // Weights
                    else if (line.find("_Weights") != std::string::npos) {
                        int layerIdx = getLayerIdx(line);
                        if (layerIdx < 0 || layerIdx >= (int)loadedNetwork.size()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Invalid layer index: " << layerIdx;
                        }
                        if (loadedNetwork[layerIdx].W.empty() || loadedNetwork[layerIdx].W[0].empty()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Weights encountered before layer size definition";
                        }
                        for (auto& row : loadedNetwork[layerIdx].W) {
                            if (!std::getline(in, line)) {
                                throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF while reading weights";
                            }
                            auto vals = split(line);
                            if (vals.size() != row.size()) {
                                throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                                    << "Weight row size mismatch for layer " << layerIdx;
                            }
                            for (size_t j = 0; j < vals.size(); ++j)
                                row[j] = std::stod(vals[j]);
                        }
                    }
                    // Biases
                    else if (line.find("_Biases") != std::string::npos) {
                        int layerIdx = getLayerIdx(line);
                        if (layerIdx < 0 || layerIdx >= (int)loadedNetwork.size()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Invalid layer index: " << layerIdx;
                        }
                        if (!std::getline(in, line)) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF while reading biases";
                        }
                        auto vals = split(line);
                        if (vals.size() != loadedNetwork[layerIdx].b.size()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                                << "Bias size mismatch for layer " << layerIdx;
                        }
                        for (size_t j = 0; j < vals.size(); ++j)
                            loadedNetwork[layerIdx].b[j] = std::stod(vals[j]);
                    }
                }

                // DATA NORMALIZATION
                else if (section == "[DATA_NORMALIZATION]") {
                    if (line.rfind("numDimensions", 0) == 0) {
                        int normalizationDim = std::stoi(split(line)[1]);
                        dataMean.resize(normalizationDim);
                        dataStdev.resize(normalizationDim);
                        normMin.resize(normalizationDim);
                        normMax.resize(normalizationDim);
                        continue;
                    }
                    if (line == "dataMean") {
                        if (!std::getline(in, line)) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF while reading dataMean";
                        }
                        auto vals = split(line);
                        if (vals.size() != dataMean.size()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Data mean size mismatch";
                        }
                        for (size_t j = 0; j < vals.size(); ++j) {
                            dataMean[j] = std::stod(vals[j]);
                        }
                        continue;
                    }
                    if (line == "dataStdev") {
                        if (!std::getline(in, line)) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF while reading dataStdev";
                        }
                        auto vals = split(line);
                        if (vals.size() != dataStdev.size()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Data stdev size mismatch";
                        }
                        for (size_t j = 0; j < vals.size(); ++j) {
                            dataStdev[j] = std::stod(vals[j]);
                        }
                        continue;
                    }
                    if (line == "normMin") {
                        if (!std::getline(in, line)) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF while reading normMin";
                        }
                        auto vals = split(line);
                        if (vals.size() != normMin.size()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Norm min size mismatch";
                        }
                        for (size_t j = 0; j < vals.size(); ++j) {
                            normMin[j] = std::stod(vals[j]);
                        }
                        continue;
                    }
                    if (line == "normMax") {
                        if (!std::getline(in, line)) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF while reading normMax";
                        }
                        auto vals = split(line);
                        if (vals.size() != normMax.size()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Norm max size mismatch";
                        }
                        for (size_t j = 0; j < vals.size(); ++j) {
                            normMax[j] = std::stod(vals[j]);
                        }
                        continue;
                    }
                }

                // TRAINING HISTORY
                else if (section == "[TRAINING_HISTORY]") {
                    if (line == "EpochNumber,Loss") continue;

                    auto tokens = split(line);
                    if (tokens.size() == 2) {
                        if (tokens[0] == "numEpochs") {
                            int numEpochs = std::stoi(tokens[1]);
                            epochLosses.reserve(numEpochs);
                        } else if (tokens[0] == "trainingSampleSize")
                        {
                            // We can store this if needed for analysis, but it is not used in model reconstruction.
                            // If needed later, parse as size_t to avoid truncation, e.g. size_t trainingSampleSize = std::stoull(tokens[1]);
                            // and apply checked conversion only when assigning to narrower types.
                        }else {
                            epochLosses.push_back(std::stod(tokens[1]));
                        }
                    }
                }

                // OPTIMIZER STATE
                else if (section == "[OPTIMIZER_STATE]") {
                    if (line.rfind("adamStep,", 0) == 0) {
                        loadedAdamStep = std::stoi(split(line)[1]);
                        optimizerStateLoaded = true;
                        continue;
                    }
                    if (line.find("_mW") != std::string::npos && line.rfind("Layer", 0) == 0) {
                        int layerIdx = getLayerIdx(line);
                        if (layerIdx >= (int)loadedNetwork.size()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "mW layer index out of range: " << layerIdx;
                        }
                        if (layerIdx >= (int)loadedMW.size()) loadedMW.resize(layerIdx + 1);
                        auto& outRows = loadedMW[layerIdx];
                        outRows.resize(loadedNetwork[layerIdx].W.size());
                        for (auto& row : outRows) {
                            if (!std::getline(in, line)) throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF reading mW";
                            auto vals = split(line);
                            row.resize(vals.size());
                            for (size_t j = 0; j < vals.size(); ++j) row[j] = std::stod(vals[j]);
                        }
                    }
                    else if (line.find("_vW") != std::string::npos && line.rfind("Layer", 0) == 0) {
                        int layerIdx = getLayerIdx(line);
                        if (layerIdx >= (int)loadedNetwork.size()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "vW layer index out of range: " << layerIdx;
                        }
                        if (layerIdx >= (int)loadedVW.size()) loadedVW.resize(layerIdx + 1);
                        auto& outRows = loadedVW[layerIdx];
                        outRows.resize(loadedNetwork[layerIdx].W.size());
                        for (auto& row : outRows) {
                            if (!std::getline(in, line)) throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF reading vW";
                            auto vals = split(line);
                            row.resize(vals.size());
                            for (size_t j = 0; j < vals.size(); ++j) row[j] = std::stod(vals[j]);
                        }
                    }
                    else if (line.find("_mb") != std::string::npos && line.rfind("Layer", 0) == 0) {
                        int layerIdx = getLayerIdx(line);
                        if (layerIdx >= (int)loadedNetwork.size()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "mb layer index out of range: " << layerIdx;
                        }
                        if (layerIdx >= (int)loadedMb.size()) loadedMb.resize(layerIdx + 1);
                        if (!std::getline(in, line)) throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF reading mb";
                        auto vals = split(line);
                        loadedMb[layerIdx].resize(vals.size());
                        for (size_t j = 0; j < vals.size(); ++j) loadedMb[layerIdx][j] = std::stod(vals[j]);
                    }
                    else if (line.find("_vb") != std::string::npos && line.rfind("Layer", 0) == 0) {
                        int layerIdx = getLayerIdx(line);
                        if (layerIdx >= (int)loadedNetwork.size()) {
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "vb layer index out of range: " << layerIdx;
                        }
                        if (layerIdx >= (int)loadedVb.size()) loadedVb.resize(layerIdx + 1);
                        if (!std::getline(in, line)) throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF reading vb";
                        auto vals = split(line);
                        loadedVb[layerIdx].resize(vals.size());
                        for (size_t j = 0; j < vals.size(); ++j) loadedVb[layerIdx][j] = std::stod(vals[j]);
                    }
                    else if (line == "dimLossEMA") {
                        if (!std::getline(in, line)) throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF reading dimLossEMA";
                        auto vals = split(line);
                        loadedDimLossEMA.resize(vals.size());
                        for (size_t j = 0; j < vals.size(); ++j) loadedDimLossEMA[j] = std::stod(vals[j]);
                    }
                    else if (line == "dimWeights") {
                        if (!std::getline(in, line)) throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF reading dimWeights";
                        auto vals = split(line);
                        loadedDimWeights.resize(vals.size());
                        for (size_t j = 0; j < vals.size(); ++j) loadedDimWeights[j] = std::stod(vals[j]);
                    }
                }

                // EMA NETWORK
                else if (section == "[EMA_NETWORK]") {
                    if (line.rfind("numLayers", 0) == 0) {
                        int numLayers = std::stoi(split(line)[1]);
                        loadedEmaNetwork.resize(numLayers);
                        emaNetworkLoaded = true;
                        continue;
                    }
                    if (line.find("_OutSize") != std::string::npos) {
                        auto tokens = split(line);
                        int outSize = std::stoi(tokens[1]);
                        int layerIdx = getLayerIdx(tokens[0]);
                        loadedEmaNetwork[layerIdx].W.resize(outSize);
                        loadedEmaNetwork[layerIdx].b.resize(outSize);
                    }
                    else if (line.find("_InSize") != std::string::npos) {
                        auto tokens = split(line);
                        int inSize = std::stoi(tokens[1]);
                        int layerIdx = getLayerIdx(tokens[0]);
                        if (loadedEmaNetwork[layerIdx].W.empty())
                            throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "EMA InSize before OutSize for layer " << layerIdx;
                        for (auto& row : loadedEmaNetwork[layerIdx].W) row.resize(inSize);
                    }
                    else if (line.find("_Weights") != std::string::npos) {
                        int layerIdx = getLayerIdx(line);
                        for (auto& row : loadedEmaNetwork[layerIdx].W) {
                            if (!std::getline(in, line)) throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF reading EMA weights";
                            auto vals = split(line);
                            if (vals.size() != row.size()) throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "EMA weight row size mismatch layer " << layerIdx;
                            for (size_t j = 0; j < vals.size(); ++j) row[j] = std::stod(vals[j]);
                        }
                    }
                    else if (line.find("_Biases") != std::string::npos) {
                        int layerIdx = getLayerIdx(line);
                        if (!std::getline(in, line)) throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Unexpected EOF reading EMA biases";
                        auto vals = split(line);
                        if (vals.size() != loadedEmaNetwork[layerIdx].b.size()) throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "EMA bias size mismatch layer " << layerIdx;
                        for (size_t j = 0; j < vals.size(); ++j) loadedEmaNetwork[layerIdx].b[j] = std::stod(vals[j]);
                    }
                }
            }

            if (dim <= 0 || conditionDim < 0 || hidden <= 0 || layers <= 0) {
                throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Invalid model parameters in file";
            }
            for (size_t l = 0; l < loadedNetwork.size(); ++l) {
                if (loadedNetwork[l].W.empty() || loadedNetwork[l].b.empty()) {
                    throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Incomplete layer data in file";
                }
                for (const auto& row : loadedNetwork[l].W) {
                    if (row.size() != loadedNetwork[l].W[0].size()) {
                        throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Inconsistent row size in layer";
                    }
                }
            }
            for (size_t l = 1; l < loadedNetwork.size(); ++l) {
                if (loadedNetwork[l].W[0].size() != loadedNetwork[l-1].W.size()) {
                    throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Layer size mismatch between layers";
                }
            }

            std::ostringstream logMsg;
            logMsg << "Model loaded successfully from " << filename << "\n";
            if (optimizerStateLoaded) {
                logMsg << "  Optimizer state (Adam moments, step=" << loadedAdamStep << ") restored — training can be resumed.\n";
            } else {
                logMsg << "  No optimizer state found — fresh optimizer state will be used.\n";
            }
            logMsg << "  DimWeightController: " << (useDimWeightController ? "enabled" : "disabled");
            if (useDimWeightController) logMsg << " (EMADecay=" << dimWeightEMADecay << ")";
            logMsg << "\n";
            logMsg << "  EMANetwork: " << (useEMANetwork ? "enabled" : "disabled");
            if (useEMANetwork) logMsg << " (decay=" << emaNetworkDecay << ", " << (emaNetworkLoaded ? "weights from checkpoint" : "initialized from network") << ")";
            logMsg << "\n";
            mf::LogInfo("ScoreBasedDiffusionModel::loadModel") << logMsg.str();

            // Reconstruct model without random weight initialization
            ScoreBasedDiffusionModel model(
                randFlat,
                randGaussQ,
                dim,
                conditionDim,
                timeEmbeddingDim,
                inputEmbeddingDim,
                conditionEmbeddingDim,
                hidden,
                layers,
                optimizerType,
                adamBeta1,
                adamBeta2,
                adamEps,
                scheduleType,
                betaMin,
                betaMax,
                cosineOffset,
                logSigMin,
                logSigMax,
                epsPrediction,
                lossWeightPower,
                batchSize,
                gradientClipThreshold,
                learningRate,
                useDimWeightController,
                dimWeightEMADecay,
                useEMANetwork,
                emaNetworkDecay,
                diffusionSteps,
                false // initializeRandomWeights
            );

            // Validate loaded network shape and dimensions match the initialized model before overwriting.
            if (loadedNetwork.size() != model.network_.size()) {
                throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                    << "Layer count mismatch: loaded " << loadedNetwork.size()
                    << " layers, expected " << model.network_.size()
                    << ". The model file may be missing or malformed in [NETWORK_PARAMETERS] (e.g. numLayers).";
            }
            for (size_t l = 0; l < loadedNetwork.size(); ++l) {
                if (loadedNetwork[l].W.empty() || loadedNetwork[l].W[0].empty() || loadedNetwork[l].b.empty()) {
                    throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                        << "Loaded layer " << l
                        << " is not fully sized before overwrite (empty weights/biases).";
                }
                // Check that loaded layer dimensions match the model's initialized dimensions
                if (loadedNetwork[l].W.size() != model.network_[l].W.size() ||
                    loadedNetwork[l].W[0].size() != model.network_[l].W[0].size()) {
                    throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                        << "Layer " << l << " weight dimension mismatch: loaded W["
                        << loadedNetwork[l].W.size() << "]["
                        << loadedNetwork[l].W[0].size() << "], expected W["
                        << model.network_[l].W.size() << "]["
                        << model.network_[l].W[0].size() << "].";
                }
                if (loadedNetwork[l].b.size() != model.network_[l].b.size()) {
                    throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                        << "Layer " << l << " bias dimension mismatch: loaded b["
                        << loadedNetwork[l].b.size() << "], expected b["
                        << model.network_[l].b.size() << "].";
                }
            }

            if ((int)dataMean.size() != (dim + conditionDim)) {
                throw cet::exception("ScoreBasedDiffusionModel::loadModel") << "Normalization parameter size mismatch";
            } // dataStdev, normMin, normMax have same dimensions

            // Restore network weights
            for (size_t l = 0; l < model.network_.size(); ++l) {
                model.network_[l].W = loadedNetwork[l].W;
                model.network_[l].b = loadedNetwork[l].b;
            }
            model.dataMean_ = dataMean;
            model.dataStdev_ = dataStdev;
            model.normMin_ = normMin;
            model.normMax_ = normMax;
            model.epochLosses_ = epochLosses;

            // Restore Adam optimizer state if present, enabling seamless training resumption
            if (optimizerStateLoaded) {
                bool momentSizeOk = (loadedMW.size() == model.network_.size() &&
                                    loadedVW.size() == model.network_.size() &&
                                    loadedMb.size() == model.network_.size() &&
                                    loadedVb.size() == model.network_.size());
                if (momentSizeOk) {
                    for (size_t l = 0; l < model.network_.size() && momentSizeOk; ++l) {
                        if (loadedMW[l].size() != model.network_[l].mW.size() ||
                            loadedVW[l].size() != model.network_[l].vW.size() ||
                            loadedMb[l].size() != model.network_[l].mb.size() ||
                            loadedVb[l].size() != model.network_[l].vb.size()) {
                            momentSizeOk = false;
                            break;
                        }
                        for (size_t i = 0; i < loadedMW[l].size() && momentSizeOk; ++i) {
                            if (loadedMW[l][i].size() != model.network_[l].mW[i].size() ||
                                loadedVW[l][i].size() != model.network_[l].vW[i].size()) {
                                momentSizeOk = false;
                            }
                        }
                    }
                }
                if (momentSizeOk) {
                    model.adamStep_ = loadedAdamStep;
                    for (size_t l = 0; l < model.network_.size(); ++l) {
                        model.network_[l].mW = loadedMW[l];
                        model.network_[l].vW = loadedVW[l];
                        model.network_[l].mb = loadedMb[l];
                        model.network_[l].vb = loadedVb[l];
                    }
                } else {
                    mf::LogWarning("ScoreBasedDiffusionModel::loadModel")
                        << "Optimizer state dimensions do not match the reconstructed model — "
                        << "falling back to fresh optimizer state. Training can continue but will not resume from the saved checkpoint.";
                }

                // Restore dim weight controller state
                if (loadedDimLossEMA.size() == static_cast<size_t>(dim) &&
                    loadedDimWeights.size() == static_cast<size_t>(dim)) {
                    model.dimLossEMA_ = loadedDimLossEMA;
                    model.dimWeights_ = loadedDimWeights;
                } else if (!loadedDimLossEMA.empty()) {
                    mf::LogWarning("ScoreBasedDiffusionModel::loadModel")
                        << "Dim weight controller state size mismatch — resetting to defaults (all-zeros EMA, all-ones weights).";
                }
            }

            // Restore EMA network weights
            if (emaNetworkLoaded && useEMANetwork) {
                bool emaOk = (loadedEmaNetwork.size() == model.emaNetwork_.size());
                for (size_t l = 0; l < loadedEmaNetwork.size() && emaOk; ++l) {
                    if (loadedEmaNetwork[l].W.size() != model.emaNetwork_[l].W.size() ||
                        (!loadedEmaNetwork[l].W.empty() && loadedEmaNetwork[l].W[0].size() != model.emaNetwork_[l].W[0].size()) ||
                        loadedEmaNetwork[l].b.size() != model.emaNetwork_[l].b.size())
                        emaOk = false;
                }
                if (emaOk) {
                    for (size_t l = 0; l < model.emaNetwork_.size(); ++l) {
                        model.emaNetwork_[l].W = loadedEmaNetwork[l].W;
                        model.emaNetwork_[l].b = loadedEmaNetwork[l].b;
                    }
                    mf::LogInfo("ScoreBasedDiffusionModel::loadModel") << "EMA network weights restored from checkpoint.";
                } else {
                    mf::LogWarning("ScoreBasedDiffusionModel::loadModel")
                        << "EMA network dimensions mismatch — re-initializing EMA from loaded network weights.";
                    for (size_t l = 0; l < model.emaNetwork_.size(); ++l) {
                        model.emaNetwork_[l].W = model.network_[l].W;
                        model.emaNetwork_[l].b = model.network_[l].b;
                    }
                }
            } else if (useEMANetwork) {
                // Old checkpoint without [EMA_NETWORK] — seed EMA from the loaded network weights
                for (size_t l = 0; l < model.emaNetwork_.size(); ++l) {
                    model.emaNetwork_[l].W = model.network_[l].W;
                    model.emaNetwork_[l].b = model.network_[l].b;
                }
                mf::LogInfo("ScoreBasedDiffusionModel::loadModel")
                    << "No [EMA_NETWORK] section found — EMA network initialized from loaded network weights.";
            }

            return model;
        } else {
            throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                << "Unrecognized file extension in \"" << filename
                << "\"; expected \".bin\" or \".csv\".";
        }
    }

    SBDMGeneratedSample ScoreBasedDiffusionModel::generateSample(
        const std::vector<double>& condition,
        bool useEMANetworkIfAvailable,
        bool useHeun,
        bool useSDE,
        int diffusionSteps,
        double sdeToOdeSigmaThreshold
    )
    {
        if (condition.size() != static_cast<size_t>(conditionDim_)) {
            throw cet::exception("ScoreBasedDiffusionModel::generateSample")
                << "Conditioning dimension mismatch: got " << condition.size()
                << ", expected " << conditionDim_;
        }

        // Use provided diffusionSteps if positive, otherwise use the model's configured value
        int steps = (diffusionSteps > 0) ? diffusionSteps : diffusionSteps_;

        // Start from pure noise
        std::vector<double> x(dim_);

        for (int i = 0; i < dim_; ++i) {
            x[i] = randGaussQ_.fire();
        }

        return reverseDiffuseFrom(std::move(x), condition, steps,
                                  useEMANetworkIfAvailable, useHeun, useSDE,
                                  steps, sdeToOdeSigmaThreshold);
    }

    SBDMGeneratedSample ScoreBasedDiffusionModel::reverseDiffuseFrom(
        std::vector<double> x,
        const std::vector<double>& condition,
        int stepStart,
        bool useEMANetworkIfAvailable,
        bool useHeun,
        bool useSDE,
        int steps,
        double sdeToOdeSigmaThreshold
    )
    {
        double sigma_safe = 1e-12; // small constant to prevent division by zero in case of very small sigma

        // Reverse diffusion process
        for (int step = stepStart - 1; step >= 0; --step) {

            double t = ((double)step + 1.0)/steps;
            double dt = 1.0/steps;
            double beta_val = beta(t); // as long as diffusionSteps_ is not too large, s should not become too small to cause numerical issues.
            bool effectiveSDE = useSDE &&
                (sdeToOdeSigmaThreshold < 0.0 || sigma(t) >= sdeToOdeSigmaThreshold);

            const auto& inferNet = (useEMANetworkIfAvailable && useEMANetwork_) ? emaNetwork_ : network_;
            if (!useHeun) {
                // Euler method (1st order)
                auto input = buildNetworkInput(x, condition, t);

                auto score = forwardInference(input, inferNet);
                if (epsPrediction_) {
                    double s = std::max(sigma(t), sigma_safe);
                    for (int i = 0; i < dim_; ++i) score[i] = -score[i] / s;
                }

                for (int i = 0; i < dim_; ++i) {
                    if (effectiveSDE) {
                        // SDE solver:
                        double drift = 0.5 * beta_val * x[i] + beta_val * score[i];
                        double noise = std::sqrt(beta_val * dt) * randGaussQ_.fire();
                        x[i] += drift * dt + noise;
                    } else {
                        // ODE solver:
                        double drift = 0.5 * beta_val * x[i] + 0.5 * beta_val * score[i];
                        x[i] += drift * dt;
                    }
                }
            } else {
                // Heun's method (2nd order) Only ODE solver, no noise added

                // sahred noise vector
                std::vector<double> dw(dim_);
                double noiseScale = std::sqrt(beta_val * dt);
                for (int i = 0; i < dim_; ++i) {
                    dw[i] = noiseScale * randGaussQ_.fire();
                }

                // k1 = f(x,t)
                auto input = buildNetworkInput(x, condition, t);
                auto score_k1 = forwardInference(input, inferNet);
                if (epsPrediction_) {
                    double s = std::max(sigma(t), sigma_safe);
                    for (int i = 0; i < dim_; ++i) score_k1[i] = -score_k1[i] / s;
                }

                std::vector<double> k1(dim_);
                for (int i = 0; i < dim_; ++i) {
                    if (effectiveSDE) {
                        // SDE solver:
                        k1[i] = 0.5 * beta_val * x[i] + beta_val * score_k1[i];
                    } else {
                        // ODE solver:
                        k1[i] = 0.5 * beta_val * x[i] + 0.5 * beta_val * score_k1[i];
                    }
                }

                // predictor
                std::vector<double> x_pred(dim_);
                for (int i = 0; i < dim_; ++i) {
                    x_pred[i] = x[i] + k1[i] * dt;
                    if (effectiveSDE) {
                        x_pred[i] += dw[i]; // add noise
                    }
                }

                // next time
                double t_next = std::max(0.0, t - dt);
                double b_next = beta(t_next);

                // k2 = f(x_pred,t_next)
                auto input_next = buildNetworkInput(x_pred, condition, t_next);
                auto score_k2 = forwardInference(input_next, inferNet);
                if (epsPrediction_) {
                    double s_next = std::max(sigma(t_next), sigma_safe);
                    for (int i = 0; i < dim_; ++i) score_k2[i] = -score_k2[i] / s_next;
                }

                std::vector<double> k2(dim_);
                for (int i = 0; i < dim_; ++i) {
                    if (effectiveSDE) {
                        // SDE solver:
                        k2[i] = 0.5 * b_next * x_pred[i] + b_next * score_k2[i];
                    } else {
                        // ODE solver:
                        k2[i] = 0.5 * b_next * x_pred[i] + 0.5 * b_next * score_k2[i];
                    }
                }

                // trapezoidal update
                for (int i = 0; i < dim_; ++i) {
                    x[i] += 0.5 * (k1[i] + k2[i]) * dt;
                    if (effectiveSDE) {
                        x[i] += dw[i]; // add noise
                    }
                }
            }
        }

        SBDMGeneratedSample generatedSample;
        generatedSample.zscore = x;
        generatedSample.value.resize(dim_);
        for (int i = 0; i < dim_; ++i) {
            generatedSample.value[i] = x[i]*dataStdev_[i] + dataMean_[i];
            // note the order in the stdev and mean is always x then cond
        }
        return generatedSample;
    }

    // Partial-reverse diagnostic: noise a normalized data sample at t0 (snapped to the
    // sampler's time grid so the first network evaluation matches the noising time), then
    // integrate the same reverse process generateSample() uses from t0 down to 0.
    SBDMGeneratedSample ScoreBasedDiffusionModel::partialReverseSample(
        const std::vector<double>& xNorm,
        const std::vector<double>& condition,
        double t0,
        bool useEMANetworkIfAvailable,
        bool useHeun,
        bool useSDE,
        int diffusionSteps,
        double sdeToOdeSigmaThreshold
    )
    {
        if (xNorm.size() != static_cast<size_t>(dim_)) {
            throw cet::exception("ScoreBasedDiffusionModel::partialReverseSample")
                << "State dimension mismatch: got " << xNorm.size() << ", expected " << dim_;
        }
        if (condition.size() != static_cast<size_t>(conditionDim_)) {
            throw cet::exception("ScoreBasedDiffusionModel::partialReverseSample")
                << "Conditioning dimension mismatch: got " << condition.size()
                << ", expected " << conditionDim_;
        }
        if (t0 <= 0.0 || t0 > 1.0) {
            throw cet::exception("ScoreBasedDiffusionModel::partialReverseSample")
                << "Invalid diffusion time t0=" << t0 << ": must be in (0,1]";
        }

        // Use provided diffusionSteps if positive, otherwise use the model's configured value
        int steps = (diffusionSteps > 0) ? diffusionSteps : diffusionSteps_;

        // Snap t0 to the sampler's grid t = stepStart/steps; the reverse loop's first
        // evaluation is at exactly that time, so the noised state sits on the grid.
        int stepStart = std::min(steps, std::max(1, static_cast<int>(std::lround(t0 * steps))));
        double tGrid = static_cast<double>(stepStart) / steps;

        std::vector<double> eps;
        auto x = addNoise(xNorm, tGrid, eps);

        return reverseDiffuseFrom(std::move(x), condition, stepStart,
                                  useEMANetworkIfAvailable, useHeun, useSDE,
                                  steps, sdeToOdeSigmaThreshold);
    }

    // One-step denoising diagnostic: noise a normalized sample at fixed t, run one forward
    // pass, and reconstruct x0_hat via Tweedie's formula for the VP forward process
    //   x_t = sqrt(alphabar(t)) * x0 + sigma(t) * eps  =>  x0_hat = (x_t - sigma * eps_hat) / sqrt(alphabar)
    SBDMGeneratedSample ScoreBasedDiffusionModel::denoiseOneStep(
        const std::vector<double>& xNorm,
        const std::vector<double>& condition,
        double t,
        bool useEMANetworkIfAvailable
    )
    {
        if (xNorm.size() != static_cast<size_t>(dim_)) {
            throw cet::exception("ScoreBasedDiffusionModel::denoiseOneStep")
                << "State dimension mismatch: got " << xNorm.size() << ", expected " << dim_;
        }
        if (condition.size() != static_cast<size_t>(conditionDim_)) {
            throw cet::exception("ScoreBasedDiffusionModel::denoiseOneStep")
                << "Conditioning dimension mismatch: got " << condition.size()
                << ", expected " << conditionDim_;
        }
        if (t <= 0.0 || t >= 1.0) {
            throw cet::exception("ScoreBasedDiffusionModel::denoiseOneStep")
                << "Invalid diffusion time t=" << t << ": must be in (0,1)";
        }

        const double eps_safe = 1e-12;

        std::vector<double> eps;
        auto xt = addNoise(xNorm, t, eps);
        double s  = std::max(sigma(t), eps_safe);
        double ab = std::sqrt(std::max(alphabar(t), eps_safe));

        auto input = buildNetworkInput(xt, condition, t);
        const auto& inferNet = (useEMANetworkIfAvailable && useEMANetwork_) ? emaNetwork_ : network_;
        auto out = forwardInference(input, inferNet);

        // Recover eps_hat from the network output: the network predicts either eps directly
        // or the score = -eps/sigma.
        std::vector<double> x0hat(dim_);
        for (int i = 0; i < dim_; ++i) {
            double epsHat = epsPrediction_ ? out[i] : -out[i] * s;
            x0hat[i] = (xt[i] - s * epsHat) / ab;
        }

        SBDMGeneratedSample result;
        result.zscore = x0hat;
        result.value.resize(dim_);
        for (int i = 0; i < dim_; ++i) {
            result.value[i] = x0hat[i] * dataStdev_[i] + dataMean_[i];
        }
        return result;
    }

} // namespace mu2e
