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

// stdlib includes
#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>

#include "Offline/MachineLearningTools/inc/ScoreBasedDiffusionModel.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/STMMC/inc/VDResamplerTransforms.hh"

typedef unsigned long VolumeId_type;

namespace mu2e {
  class VDResamplerTrain : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> StepPointMCsTag{   Name("StepPointMCsTag"),         Comment("Tag identifying the StepPointMCs")};
        fhicl::Atom<art::InputTag> SimParticlemvTag{  Name("SimParticlemvTag"),        Comment("Tag identifying the SimParticlemv")};
        fhicl::Atom<bool> SBDMuseTwoStageTraining{    Name("SBDMuseTwoStageTraining"), Comment("If true, train the two-stage factorized model. If false, train all 6 dimensions at once."), true};
        fhicl::Atom<std::string> SBDMallAtOnceModelFile{ Name("SBDMallAtOnceModelFile"), Comment("CSV filename for the all-at-once 6D SBDM model parameters"), ""};
        fhicl::Atom<std::string> SBDMstage1ModelFile{ Name("SBDMstage1ModelFile"),     Comment("CSV filename for the trained stage-1 SBDM model parameters"), ""};
        fhicl::Atom<std::string> SBDMstage2ModelFile{ Name("SBDMstage2ModelFile"),     Comment("CSV filename for the trained stage-2 SBDM model parameters"), ""};
        fhicl::Atom<int> VirtualDetectorID{           Name("VirtualDetectorID"),       Comment("ID of the virtual detector to train on"),                        116};
        fhicl::Atom<double> VDz0{                     Name("VDz0"),                    Comment("z coordinate of the virtual detector"),                          37700.39};
        fhicl::Atom<double> VDr{                      Name("VDr"),                     Comment("VD radius"),                                                     2000.0 };
        fhicl::Atom<int> pdgID{                       Name("pdgID"),                   Comment("pdgID of the particle to train on"),                             22};
        fhicl::Atom<int> SBDMhidden{                  Name("SBDMhidden"),              Comment("Size of hidden layers in the SBDM neural network"),              128};
        fhicl::Atom<int> SBDMlayers{                  Name("SBDMlayers"),              Comment("Number of layers in the SBDM neural network"),                   4};
        fhicl::Atom<std::string> SBDMoptimizer{       Name("SBDMoptimizer"),           Comment("Optimizer for training the SBDM neural network (SGD or ADAM)"), "ADAM"};
        fhicl::Atom<double> SBDMadamBeta1{            Name("SBDMadamBeta1"),           Comment("Adam optimizer beta1 parameter"),                                0.9};
        fhicl::Atom<double> SBDMadamBeta2{            Name("SBDMadamBeta2"),           Comment("Adam optimizer beta2 parameter"),                                0.999};
        fhicl::Atom<double> SBDMadamEps{              Name("SBDMadamEps"),             Comment("Adam optimizer epsilon parameter"),                              1e-8};
        fhicl::Atom<std::string> SBDMnoiseSchedule{   Name("SBDMnoiseSchedule"),       Comment("Noise schedule for the SBDM (LINEAR or COSINE)"),                "COSINE"};
        fhicl::Atom<double> SBDMbetaMin{              Name("SBDMbetaMin"),             Comment("Minimum noise schedule parameter (for LINEAR schedule)") ,       1e-4};
        fhicl::Atom<double> SBDMbetaMax{              Name("SBDMbetaMax"),             Comment("Maximum noise schedule parameter (for LINEAR schedule)"),        0.02};
        fhicl::Atom<double> SBDMcosineOffset{         Name("SBDMcosineOffset"),        Comment("Offset parameter (for cosine schedule)"),                        0.008};
        fhicl::Atom<int> SBDMbatchSize{               Name("SBDMbatchSize"),           Comment("Batch size for training the SBDM"),                              32};
        fhicl::Atom<double> SBDMgradientClip{         Name("SBDMgradientClip"),        Comment("Gradient clipping threshold for training the SBDM"),             1.0};
        fhicl::Atom<double> SBDMlearningRate{         Name("SBDMlearningRate"),        Comment("Learning rate for training the SBDM"),                           1e-3};
        fhicl::Atom<int> SBDMdiffusionSteps{          Name("SBDMdiffusionSteps"),      Comment("Number of steps in the diffusion process for the SBDM"),         200};
        fhicl::Atom<int> SBDMtrainingSize{            Name("SBDMtrainingSize"),        Comment("Size of the training data for the SBDM"),                        -1}; // -1 means use all available data
        fhicl::Atom<int> SBDMtrainingEpochs{          Name("SBDMtrainingEpochs"),      Comment("Number of epochs to train the selected SBDM mode"),              10};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit VDResamplerTrain(const Parameters& conf);
      void analyze(const art::Event& e);
      void endJob();
    private:
      art::RandomNumberGenerator::base_engine_t& engine;
      CLHEP::RandFlat randFlat_;
      CLHEP::RandGaussQ randGaussQ_;
      art::ProductToken<StepPointMCCollection> StepPointMCsToken;
      art::ProductToken<SimParticleCollection> SimParticlemvToken;

      // SBDM model + data
      bool useTwoStageTraining = true;
      std::unique_ptr<ScoreBasedDiffusionModel> allAtOnceModel;
      std::unique_ptr<ScoreBasedDiffusionModel> stage1Model;
      std::unique_ptr<ScoreBasedDiffusionModel> stage2Model;
      std::vector<DiffusionTrainingSample> allAtOnceTrainingData;
      std::vector<DiffusionTrainingSample> stage1TrainingData;
      std::vector<DiffusionTrainingSample> stage2TrainingData;
      std::string SBDMallAtOnceModelFile;
      std::string SBDMstage1ModelFile;
      std::string SBDMstage2ModelFile;
      VolumeId_type VirtualDetectorID = 0;
      double VDz0 = 0.0;
      double VDr = 0.0;
      int pdgID = 0;
      int trainingEpochs = 0;
      int trainingSize = -1;

      // variables to read from the art event
      int stepPdgId = 0;
      double x = 0.0, y = 0.0, z = 0.0, time = 0.0;
      double px = 0.0, py = 0.0, pz = 0.0;
      VolumeId_type virtualdetectorId = 0;

      // transform variables for training data preparation
      double x0 = VDResampler::kX0;
      double y0 = VDResampler::kY0;
      // time scaling
      double t0 = VDResampler::kT0;
      double tScale = VDResampler::kTScale;
      // momentum scaling
      double p0 = VDResampler::kP0;
  };

  VDResamplerTrain::VDResamplerTrain(const Parameters& conf) :
    art::EDAnalyzer(conf),
    engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    randFlat_(engine),
    randGaussQ_(engine),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag())),
    useTwoStageTraining(conf().SBDMuseTwoStageTraining()),
    SBDMallAtOnceModelFile(conf().SBDMallAtOnceModelFile()),
    SBDMstage1ModelFile(conf().SBDMstage1ModelFile()),
    SBDMstage2ModelFile(conf().SBDMstage2ModelFile()),
    VirtualDetectorID(conf().VirtualDetectorID()),
    VDz0(conf().VDz0()),
    VDr(conf().VDr()),
    pdgID(conf().pdgID()),
    trainingEpochs(conf().SBDMtrainingEpochs()),
    trainingSize(conf().SBDMtrainingSize())
  {
    // Validate geometry configuration
    if (VDr <= 0.0) {
        throw cet::exception("VDResamplerTrain")
            << "VDr must be positive (got " << VDr << "); "
            << "rho = r/VDr would produce inf/NaN in training data.";
    }
    if (!std::isfinite(VDz0)) {
        throw cet::exception("VDResamplerTrain")
            << "VDz0 must be finite (got " << VDz0 << ").";
    }

    // optimizer selection
    ScoreBasedDiffusionModel::OptimizerType opt;
    if (conf().SBDMoptimizer() == "SGD") {
        opt = ScoreBasedDiffusionModel::OptimizerType::SGD;
    } else {
        if (conf().SBDMoptimizer() != "ADAM") {
            mf::LogWarning("VDResamplerTrain")
                << "Unrecognized SBDMoptimizer value \"" << conf().SBDMoptimizer()
                << "\"; falling back to ADAM.";
        }
        opt = ScoreBasedDiffusionModel::OptimizerType::ADAM;
    }

    // noise schedule
    ScoreBasedDiffusionModel::NoiseScheduleType sched;
    if (conf().SBDMnoiseSchedule() == "LINEAR") {
        sched = ScoreBasedDiffusionModel::NoiseScheduleType::LINEAR;
    } else {
        if (conf().SBDMnoiseSchedule() != "COSINE") {
            mf::LogWarning("VDResamplerTrain")
                << "Unrecognized SBDMnoiseSchedule value \"" << conf().SBDMnoiseSchedule()
                << "\"; falling back to COSINE.";
        }
        sched = ScoreBasedDiffusionModel::NoiseScheduleType::COSINE;
    }

    if (useTwoStageTraining) {
      // create stage-1 model for (t', x', y')
      stage1Model = std::make_unique<ScoreBasedDiffusionModel>(
          randFlat_,
          randGaussQ_,
          3, // dim
          0, // conditionDim
          conf().SBDMhidden(),
          conf().SBDMlayers(),
          opt,
          conf().SBDMadamBeta1(),
          conf().SBDMadamBeta2(),
          conf().SBDMadamEps(),
          sched,
          conf().SBDMbetaMin(),
          conf().SBDMbetaMax(),
          conf().SBDMcosineOffset(),
          conf().SBDMbatchSize(),
          conf().SBDMgradientClip(),
          conf().SBDMlearningRate(),
          conf().SBDMdiffusionSteps()
      );

      // create stage-2 model for (p_r', p_phi', p_z' | t', x', y')
      stage2Model = std::make_unique<ScoreBasedDiffusionModel>(
          randFlat_,
          randGaussQ_,
          3, // dim
          3, // conditionDim
          conf().SBDMhidden(),
          conf().SBDMlayers(),
          opt,
          conf().SBDMadamBeta1(),
          conf().SBDMadamBeta2(),
          conf().SBDMadamEps(),
          sched,
          conf().SBDMbetaMin(),
          conf().SBDMbetaMax(),
          conf().SBDMcosineOffset(),
          conf().SBDMbatchSize(),
          conf().SBDMgradientClip(),
          conf().SBDMlearningRate(),
          conf().SBDMdiffusionSteps()
      );
      // allocate memory according to the training size (if specified, otherwise will grow dynamically)
      if (trainingSize > 0) {
        stage1TrainingData.reserve(trainingSize);
        stage2TrainingData.reserve(trainingSize);
      } else {
        stage1TrainingData.reserve(1000);
        stage2TrainingData.reserve(1000);
      }
    } else {
      allAtOnceModel = std::make_unique<ScoreBasedDiffusionModel>(
          randFlat_,
          randGaussQ_,
          6, // dim
          0, // conditionDim
          conf().SBDMhidden(),
          conf().SBDMlayers(),
          opt,
          conf().SBDMadamBeta1(),
          conf().SBDMadamBeta2(),
          conf().SBDMadamEps(),
          sched,
          conf().SBDMbetaMin(),
          conf().SBDMbetaMax(),
          conf().SBDMcosineOffset(),
          conf().SBDMbatchSize(),
          conf().SBDMgradientClip(),
          conf().SBDMlearningRate(),
          conf().SBDMdiffusionSteps()
      );
      // allocate memory according to the training size (if specified, otherwise will grow dynamically)
      if (trainingSize > 0) {
        allAtOnceTrainingData.reserve(trainingSize);
      } else {
        allAtOnceTrainingData.reserve(1000);
      }
    }
  };

  void VDResamplerTrain::analyze(const art::Event& event) {
    // Get the data products from the event
    auto const& StepPointMCs = event.getProduct(StepPointMCsToken);
    if (StepPointMCs.empty())
      return;
    auto const& SimParticles = event.getProduct(SimParticlemvToken);
    if (SimParticles.empty())
      return;

    // Loop over all VD hits
    for (const StepPointMC& step : StepPointMCs) {
      // Get the associated particle
      const SimParticle& particle = SimParticles.at(step.trackId());

      // Extract the parameters
      stepPdgId = particle.pdgId();
      virtualdetectorId = step.virtualDetectorId();
      time = step.time();
      x = step.position().x(); // This coordinate is in Mu2e frame, will be shifted to relative x w.r.t. the beamline
      y = step.position().y();
      z = step.position().z();
      px = step.momentum().x();
      py = step.momentum().y();
      pz = step.momentum().z();

      if (virtualdetectorId != VirtualDetectorID || (stepPdgId != pdgID && pdgID != 0) || pz <= 0)
        continue; // Filter hits based on the virtual detector ID, particle type, and pz

      double x_trans = 0.0;
      double y_trans = 0.0;
      double t_trans = 0.0;
      double pr_t = 0.0;
      double pphi_t = 0.0;
      double pz_t = 0.0;
      VDResampler::forwardTransformSample(
        x,
        y,
        z,
        time,
        px,
        py,
        pz,
        x0,
        y0,
        t0,
        tScale,
        p0,
        VDr,
        VDz0,
        x_trans,
        y_trans,
        t_trans,
        pr_t,
        pphi_t,
        pz_t
      );

      if (useTwoStageTraining) {
        DiffusionTrainingSample stage1Sample;
        stage1Sample.x = {t_trans, x_trans, y_trans};
        stage1TrainingData.push_back(std::move(stage1Sample));

        DiffusionTrainingSample stage2Sample;
        stage2Sample.x = {pr_t, pphi_t, pz_t};
        stage2Sample.cond = {t_trans, x_trans, y_trans};
        stage2TrainingData.push_back(std::move(stage2Sample));
      } else {
        DiffusionTrainingSample allAtOnceSample;
        allAtOnceSample.x = {t_trans, x_trans, y_trans, pr_t, pphi_t, pz_t};
        allAtOnceTrainingData.push_back(std::move(allAtOnceSample));
      }
    };
    return;
  };

  void VDResamplerTrain::endJob() {

      if (useTwoStageTraining && (stage1TrainingData.empty() || stage2TrainingData.empty())) {
        mf::LogWarning("VDResamplerTrain") << "No training data collected.";
        return;
      }
      if (!useTwoStageTraining && allAtOnceTrainingData.empty()) {
        mf::LogWarning("VDResamplerTrain") << "No training data collected.";
        return;
      }

      // if SBDMtrainingSize is set and smaller than the collected training data,
      // truncate the training data to the specified size.
      if (useTwoStageTraining) {
        if(trainingSize > 0 && (int)stage1TrainingData.size() > trainingSize)
          stage1TrainingData.resize(trainingSize);
        if(trainingSize > 0 && (int)stage2TrainingData.size() > trainingSize)
          stage2TrainingData.resize(trainingSize);

        if (SBDMstage1ModelFile.empty() || SBDMstage2ModelFile.empty()) {
          throw cet::exception("VDResamplerTrain") << "Two-stage training requires both SBDMstage1ModelFile and SBDMstage2ModelFile.";
        }

        mf::LogInfo("VDResamplerTrain")
            << "Training stage-1 diffusion model with " << stage1TrainingData.size()
            << " samples and " << trainingEpochs << " epochs...";
        stage1Model->train(stage1TrainingData, trainingEpochs);
        stage1Model->saveModel(SBDMstage1ModelFile);
        mf::LogInfo("VDResamplerTrain") << "Stage-1 model saved to " << SBDMstage1ModelFile;

        mf::LogInfo("VDResamplerTrain")
            << "Training stage-2 diffusion model with " << stage2TrainingData.size()
            << " samples and " << trainingEpochs << " epochs...";
        stage2Model->train(stage2TrainingData, trainingEpochs);
        stage2Model->saveModel(SBDMstage2ModelFile);
        mf::LogInfo("VDResamplerTrain") << "Stage-2 model saved to " << SBDMstage2ModelFile;

      } else {
        if(trainingSize > 0 && (int)allAtOnceTrainingData.size() > trainingSize)
          allAtOnceTrainingData.resize(trainingSize);

        if (SBDMallAtOnceModelFile.empty()) {
          throw cet::exception("VDResamplerTrain") << "All-at-once training requires SBDMallAtOnceModelFile.";
        }

        mf::LogInfo("VDResamplerTrain")
            << "Training all-at-once diffusion model with " << allAtOnceTrainingData.size()
            << " samples and " << trainingEpochs << " epochs...";
        allAtOnceModel->train(allAtOnceTrainingData, trainingEpochs);
        allAtOnceModel->saveModel(SBDMallAtOnceModelFile);

        mf::LogInfo("VDResamplerTrain") << "All-at-once model saved to " << SBDMallAtOnceModelFile;
      }

     return;
  };

}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::VDResamplerTrain)
