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

    // Linear interpolation of beta(t) between betaMin_ and betaMax_.
    // This is only used in the linear noise schedule.
    double ScoreBasedDiffusionModel::beta(double t) const {
        return betaMin_ + t * (betaMax_ - betaMin_);
    }

    // Standard deviation of noise at diffusion time t.
    // For the linear schedule, use simply sqrt(beta(t)).
    // For the cosine schedule, this is derived from the cumulative noise schedule.
    double ScoreBasedDiffusionModel::sigma(double t) const {
        if (noiseScheduleType_ == NoiseScheduleType::COSINE) {
            // Cosine noise schedule
            double f = (t + cosineOffset_) / (1.0 + cosineOffset_);
            double alpha_bar = std::cos(f * M_PI * 0.5);
            alpha_bar *= alpha_bar;
            return std::sqrt(1.0 - alpha_bar);
        } else {
            // Linear noise schedule (default)
            return std::sqrt(beta(t));
        }
    }


    ScoreBasedDiffusionModel::ScoreBasedDiffusionModel(
        CLHEP::RandFlat& randFlat,
        CLHEP::RandGaussQ& randGaussQ,
        int dim,
        int conditionDim,
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
        int batchSize,
        double gradientClipThreshold,
        double learningRate,
        int diffusionSteps
    ) : randFlat_(randFlat), randGaussQ_(randGaussQ),
        dim_(dim), conditionDim_(conditionDim), hidden_(hidden), layers_(layers),
        optimizerType_(optimizerType),
        adamBeta1_(adamBeta1), adamBeta2_(adamBeta2), adamEps_(adamEps),
        noiseScheduleType_(scheduleType),
        betaMin_(betaMin), betaMax_(betaMax), cosineOffset_(cosineOffset),
        batchSize_(batchSize), gradientClipThreshold_(gradientClipThreshold), learningRate_(learningRate),
        diffusionSteps_(diffusionSteps),
        runningLoss_(0.0), adamStep_(0), trainingSampleSize_(0), epochLosses_() {

        // Validate model dimensions and parameters
        if (dim <= 0 || conditionDim < 0 || hidden <= 0 || layers <= 0) {
            throw cet::exception("ScoreBasedDiffusionModel::initialization") << "Invalid model dimensions";
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
        // Input dimension = dim_ + conditionDim_ + 1
        // (+1 because diffusion time t is appended to the input vector)
        // ------------------------------------------------------------

        int inputSize = dim_ + conditionDim_ + 1;
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
            layer.W.resize(out, std::vector<double>(in));

            // Allocate biases
            layer.b.resize(out);

            // Allocate gradient buffers
            layer.gradW.resize(out, std::vector<double>(in, 0.0));
            layer.gradb.resize(out, 0.0);

            // Allocate Adam buffers
            layer.mW.resize(out, std::vector<double>(in, 0.0));
            layer.vW.resize(out, std::vector<double>(in, 0.0));

            layer.mb.resize(out, 0.0);
            layer.vb.resize(out, 0.0);

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

        // Print layer and model configuration
        std::ostringstream oss;
        oss << "ScoreBasedDiffusionModel initialized\n"
            << "Model configuration:\n"
            << "  Network architecture:\n"
            << "    | dim=" << dim_ << "\n"
            << "    | conditionDim=" << conditionDim_ << "\n"
            << "    | hidden=" << hidden_ << "\n"
            << "    | layers=" << layers_ << "\n"
            << "  Optimizer configuration:\n"
            << "    | Optimizer=" << (optimizerType_ == OptimizerType::ADAM ? "Adam" : "SGD") << "\n";
        if (optimizerType_ == OptimizerType::ADAM) {
            oss << "      |- AdamBeta1=" << adamBeta1_ << "\n"
                << "      |- AdamBeta2=" << adamBeta2_ << "\n"
                << "      |- AdamEps=" << adamEps_ << "\n";
        }
        oss << "  Noise schedule configuration:\n"
            << "    | NoiseSchedule=" << (noiseScheduleType_ == NoiseScheduleType::COSINE ? "Cosine" : "Linear") << "\n";
        if (noiseScheduleType_ == NoiseScheduleType::COSINE) {
            oss << "      |- CosineOffset=" << cosineOffset_ << "\n";
        } else {
            oss << "      |- BetaMin=" << betaMin_ << "\n"
                << "      |- BetaMax=" << betaMax_ << "\n";
        }
        oss << "  Training configuration:\n"
            << "    | BatchSize=" << batchSize_ << "\n"
            << "    | GradientClipThreshold=" << gradientClipThreshold_ << "\n"
            << "    | LearningRate=" << learningRate_ << "\n"
            << "  Diffusion process configuration:\n"
            << "    | DiffusionSteps=" << diffusionSteps_ << "\n";
        mf::LogInfo("ScoreBasedDiffusionModel::initialize") << oss.str();

        // Reserve space for forward-pass caches
        activations_.reserve(layers_ + 1);
        preactivations_.reserve(layers_);
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
            if(l != network_.size()-1) // No activation on the last layer (output layer), as it predicts the score directly.
            {
                for (size_t i=0;i<z.size();++i)
                    y[i] = silu(z[i]);
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
            if(l == static_cast<int>(network_.size())-1)
            {
                gradZ = grad;
            }
            // For hidden layers, apply the derivative of the activation function (SiLU) to the gradient.
            else
            {
                for (size_t i=0;i<grad.size();++i)
                    gradZ[i] = grad[i]*siluDeriv(z[i]);
            }

            // Compute gradients w.r.t. weights and biases for the current layer
            for (size_t i = 0; i < layer.W.size(); i++)
            {
                for (size_t j=0;j<layer.W[i].size();++j)
                {
                    // Gradient of loss w.r.t. weight W[i][j]
                    layer.gradW[i][j] += gradZ[i]*aPrev[j];
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
            std::vector<double> gradPrev(layer.W[0].size(),0.0);

            for (size_t j=0;j<gradPrev.size();++j)
                for (size_t i=0;i<layer.W.size();++i)
                    gradPrev[j] += layer.W[i][j]*gradZ[i];

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

            for (size_t i=0;i<layer.W.size();++i)
            {
                // Update weights using Adam update rule
                for (size_t j=0;j<layer.W[i].size();++j)
                {
                    double g = layer.gradW[i][j];

                    layer.mW[i][j] = adamBeta1_*layer.mW[i][j] + (1-adamBeta1_)*g;
                    layer.vW[i][j] = adamBeta2_*layer.vW[i][j] + (1-adamBeta2_)*g*g;

                    double mhat = layer.mW[i][j]/(1-std::pow(adamBeta1_,adamStep_));
                    double vhat = layer.vW[i][j]/(1-std::pow(adamBeta2_,adamStep_));

                    layer.W[i][j] -= lr*mhat/(std::sqrt(vhat)+adamEps_);

                    layer.gradW[i][j] = 0.0;
                }

                // Update bias using Adam update rule
                double g = layer.gradb[i];

                layer.mb[i] = adamBeta1_*layer.mb[i] + (1-adamBeta1_)*g;
                layer.vb[i] = adamBeta2_*layer.vb[i] + (1-adamBeta2_)*g*g;

                double mhat = layer.mb[i]/(1-std::pow(adamBeta1_,adamStep_));
                double vhat = layer.vb[i]/(1-std::pow(adamBeta2_,adamStep_));

                layer.b[i] -= lr*mhat/(std::sqrt(vhat)+adamEps_);

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

        double s = sigma(t);

        // Generate Gaussian noise vector eps of dimension dim_ using the external random engine.
        eps.resize(dim_);

        std::vector<double> xt(dim_);

        for (int i = 0; i < dim_; ++i) {
            eps[i] = randGaussQ_.fire();   // Gaussian N(0,1)
            xt[i]  = x[i] + s * eps[i];
        }

        return xt;
    }

    // Compute Mean Squared Error per dimension between predicted score and target score.
    double ScoreBasedDiffusionModel::computeLoss(
        const std::vector<double>& score,
        const std::vector<double>& target
    ) const{

        double loss = 0.0;

        for (int i=0;i<dim_;++i){
            double d = score[i]-target[i];
            loss += d*d;
        }

        return loss/dim_;
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

        // If the norm is below the threshold, no clipping is needed.
        if(norm <= maxNorm)
            return;

        // Scale down the gradients to have the specified maximum norm.
        double scale = maxNorm / norm;
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
    void ScoreBasedDiffusionModel::train(
        const std::vector<DiffusionTrainingSample>& data,
        int epochs
    )
    {
        // Check that the network is properly initialized before training.
        if (network_.empty() || network_[0].W.empty() || network_[0].W[0].empty()) {
            throw cet::exception("ScoreBasedDiffusionModel::train")
                << "Network is not properly initialized";
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
                std::swap(indices[i], indices[j]);
            }

            // Counter for number of samples processed in the epoch (used for averaging loss). Avoid using N directly to
            // allow for early stopping or partial epoch processing if needed.
            int n = 0;

            int batchCounter = 0;
            double epochLoss = 0.0;

            // iterate over the shuffled data samples
            for (size_t idx = 0; idx < N; ++idx)
            {
                const auto& sample = data[indices[idx]];

                // Sample diffusion time
                double t = randFlat_.fire();
                // container for noise vector eps
                std::vector<double> eps;

                auto xt = addNoise(sample.x,t,eps);
                double s = sigma(t);

                // The target score is the negative of the noise scaled by the noise standard deviation, i.e., -eps/sigma(t).
                std::vector<double> target(dim_);
                for (int i=0;i<dim_;++i)
                    target[i] = -eps[i]/std::max(s, eps_safe); // add eps_safe to prevent division by zero in case of very small sigma

                // Prepare input vector for the network by concatenating the noisy sample xt with the diffusion time t.
                std::vector<double> input = xt;
                input.insert(input.end(), sample.cond.begin(), sample.cond.end());
                input.push_back(t);

                // Check that input dimension matches expected dimension (dim_ + conditionDim_ + 1)
                if (input.size() != network_[0].W[0].size()) {
                    throw cet::exception("ScoreBasedDiffusionModel::train")
                        << "Training input dimension mismatch: got " << input.size()
                        << ", expected " << network_[0].W[0].size();
                }

                // Forward pass to compute the predicted score from the network given the input (noisy sample + time).
                auto score = forward(input);

                // Compute the loss (Mean Squared Error) between the predicted score and the target score.
                std::vector<double> grad(dim_);
                double loss=0.0;
                for (int i=0;i<dim_;++i)
                {
                    double d = score[i]-target[i];
                    loss += d*d;
                    // Gradient of the loss w.r.t. the predicted score is 2*(score-target)/dim_ (the division by dim_ is for averaging the loss per dimension).
                    grad[i] = 2*d/dim_;
                }

                // Backward pass to compute gradients of the loss w.r.t. network parameters using the computed gradient of the loss w.r.t. the predicted score.
                backward(grad);

                // Increment batch counter and apply optimizer step if batch size is reached. This allows for vectorized updates after processing a batch of
                // samples, which can improve training efficiency and convergence.
                batchCounter++;
                if(batchCounter == batchSize_)
                {
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
                    batchCounter = 0;
                }

                epochLoss += loss; // Accumulate loss for monitoring.
                n++;
            }

            // Final flush: apply optimizer to remaining gradients
            if(batchCounter > 0)
            {
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
            }

            epochLoss /= n;
            mf::LogInfo("ScoreBasedDiffusionModel::train") << "Epoch " << e << "  Loss=" << epochLoss;
            epochLosses_.push_back(epochLoss);
        }
    }

    void ScoreBasedDiffusionModel::saveModel(
        const std::string& filename
    )
    {
        std::ofstream out(filename);

        if (!out) {
            throw cet::exception("ScoreBasedDiffusionModel::saveModel")
                << "Cannot open file " << filename;
        }

        // Write CSV header and model configuration parameters with annotations
        out << "[MODEL_PARAMETERS]\n";
        out << "Parameter,Value\n";
        // Model architecture
        out << "dim," << dim_ << "\n";
        out << "conditionDim," << conditionDim_ << "\n";
        out << "hidden," << hidden_ << "\n";
        out << "layers," << layers_ << "\n";
        // Optimizer configuration
        out << "optimizerType," << (optimizerType_ == OptimizerType::ADAM ? "ADAM" : "SGD") << "\n";
        out << "adamBeta1," << adamBeta1_ << "\n";
        out << "adamBeta2," << adamBeta2_ << "\n";
        out << "adamEps," << adamEps_ << "\n";
        // Noise schedule configuration
        out << "noiseScheduleType," << (noiseScheduleType_ == NoiseScheduleType::COSINE ? "COSINE" : "LINEAR") << "\n";
        out << "betaMin," << betaMin_ << "\n";
        out << "betaMax," << betaMax_ << "\n";
        out << "cosineOffset," << cosineOffset_ << "\n";
        // Training configuration
        out << "batchSize," << batchSize_ << "\n";
        out << "gradientClipThreshold," << gradientClipThreshold_ << "\n";
        out << "learningRate," << learningRate_ << "\n";
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

        // Write training history
        out << "\n[TRAINING_HISTORY]\n";
        out << "numEpochs," << epochLosses_.size() << "\n";
        out << "trainingSampleSize," << trainingSampleSize_ << "\n";
        out << "EpochNumber,Loss\n";
        for (size_t i = 0; i < epochLosses_.size(); ++i) {
            out << i << "," << epochLosses_[i] << "\n";
        }
    }

    ScoreBasedDiffusionModel ScoreBasedDiffusionModel::loadModel(
        CLHEP::RandFlat& randFlat,
        CLHEP::RandGaussQ& randGaussQ,
        const std::string& filename
    )
    {
        std::ifstream in(filename);

        if (!in) {
            throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                << "Cannot open file " << filename;
        }

        std::string line;

        // Temporary storage
        int dim = 0, conditionDim = 0, hidden = 0, layers = 0;
        OptimizerType optimizerType = OptimizerType::ADAM;
        double adamBeta1 = 0.0, adamBeta2 = 0.0, adamEps = 0.0;
        NoiseScheduleType scheduleType = NoiseScheduleType::COSINE;
        double betaMin = 0.0, betaMax = 0.0, cosineOffset = 0.0;
        int batchSize = 1, diffusionSteps = 1;
        double gradientClipThreshold = 0.0, learningRate = 0.0;

        std::vector<Layer> loadedNetwork;
        std::vector<double> epochLosses;

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
                else if (key == "hidden") hidden = std::stoi(val);
                else if (key == "layers") layers = std::stoi(val);
                else if (key == "optimizerType")
                    optimizerType = (val == "ADAM") ? OptimizerType::ADAM : OptimizerType::SGD;
                else if (key == "adamBeta1") adamBeta1 = std::stod(val);
                else if (key == "adamBeta2") adamBeta2 = std::stod(val);
                else if (key == "adamEps") adamEps = std::stod(val);
                else if (key == "noiseScheduleType")
                    scheduleType = (val == "COSINE") ? NoiseScheduleType::COSINE : NoiseScheduleType::LINEAR;
                else if (key == "betaMin") betaMin = std::stod(val);
                else if (key == "betaMax") betaMax = std::stod(val);
                else if (key == "cosineOffset") cosineOffset = std::stod(val);
                else if (key == "batchSize") batchSize = std::stoi(val);
                else if (key == "gradientClipThreshold") gradientClipThreshold = std::stod(val);
                else if (key == "learningRate") learningRate = std::stod(val);
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
                    if (layerIdx == 0 && inSize != dim + conditionDim + 1) {
                        throw cet::exception("ScoreBasedDiffusionModel::loadModel")
                            << "Input layer input size mismatch (expected "
                            << (dim + conditionDim + 1) << ", got " << inSize << ")";
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
                        // We can store this if needed for analysis, but it is not used in model reconstruction, so we will ignore it for now.
                    }else {
                        epochLosses.push_back(std::stod(tokens[1]));
                    }
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
        logMsg << "Model parameters loaded successfully from " << filename << "\n"
               << "Warning: optimizer state (e.g., Adam moments) is not saved/loaded.\n"
               << "The loaded model is suitable for inference, or for retraining with a fresh optimizer state,\n"
               << "but does NOT resume training from the original state.";
        mf::LogInfo("ScoreBasedDiffusionModel::loadModel") << logMsg.str();

        // Reconstruct model
        ScoreBasedDiffusionModel model(
            randFlat,
            randGaussQ,
            dim,
            conditionDim,
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
            batchSize,
            gradientClipThreshold,
            learningRate,
            diffusionSteps
        );

        // Overwrite initialized weights with loaded values
        for (size_t l = 0; l < model.network_.size(); ++l) {
            model.network_[l].W = loadedNetwork[l].W;
            model.network_[l].b = loadedNetwork[l].b;
        }
        model.epochLosses_ = epochLosses;

        return model;
    }

    std::vector<double> ScoreBasedDiffusionModel::generateSample(
        const std::vector<double>& condition,
        bool useHeun,
        int diffusionSteps
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

        for (int i=0;i<dim_;++i)
            x[i] = randGaussQ_.fire();

        // Reverse diffusion process
        for (int step = steps - 1; step >= 0; --step) {

            double t = (double)step/steps;
            double dt = 1.0/steps;
            double s = sigma(t); // as long as diffusionSteps_ is not too large, s should not become too small to cause numerical issues.

            if (!useHeun) {
                // Euler method (1st order)
                std::vector<double> input = x;
                input.insert(input.end(), condition.begin(), condition.end());
                input.push_back(t);

                auto score = forward(input);

                for (int i=0;i<dim_;++i){
                    x[i] += -s*s*score[i]*dt;
                }
            } else {
                // Heun's method (2nd order)
                // k1 = f(x,t)
                std::vector<double> input = x;
                input.insert(input.end(), condition.begin(), condition.end());
                input.push_back(t);
                auto score_k1 = forward(input);

                std::vector<double> k1(dim_);
                for (int i=0;i<dim_;++i)
                    k1[i] = -s*s*score_k1[i];

                // predictor
                std::vector<double> x_pred(dim_);
                for (int i=0;i<dim_;++i)
                    x_pred[i] = x[i] + k1[i]*dt;

                // next time
                double t_next = std::max(0.0, t - dt);
                double s_next = sigma(t_next);

                // k2 = f(x_pred,t_next)
                std::vector<double> input_next = x_pred;
                input_next.insert(input_next.end(), condition.begin(), condition.end());
                input_next.push_back(t_next);
                auto score_k2 = forward(input_next);

                std::vector<double> k2(dim_);
                for (int i=0;i<dim_;++i)
                    k2[i] = -s_next*s_next*score_k2[i];

                // trapezoidal update
                for (int i=0;i<dim_;++i)
                    x[i] += 0.5*(k1[i]+k2[i])*dt;
            }
        }

        return x;
    }

} // namespace mu2e
