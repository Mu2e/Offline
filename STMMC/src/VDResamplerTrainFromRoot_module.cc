// VDResamplerTrainFromRoot_module.cc
// A copy of VDResamplerTrain_module.cc but takes ROOT file input instead of 
// StepPointMCCollection and SimParticleCollection from art::Event.
// Using ROOT file input from VDResamplerConfigure_module.cc,
// and saving models at specified epochs as well as the final epoch.
// Yongyi Wu, May 2026

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <string>
#include <numeric>

#include "Offline/MachineLearningTools/inc/ScoreBasedDiffusionModel.hh"
#include "Offline/STMMC/inc/VDResamplerTransforms.hh"
#include "Offline/SeedService/inc/SeedService.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "TFile.h"
#include "TTree.h"

typedef unsigned long VolumeId_type;

namespace mu2e {
  class VDResamplerTrainFromRoot : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<std::string> InputRootFile{ Name("InputRootFile"), Comment("Input ROOT file with TTree") };
        fhicl::Atom<std::string> TreeName{ Name("TreeName"), Comment("Name of the TTree in the ROOT file"), "VDResamplerTrainingSetup/ttree" };
        fhicl::Atom<bool> SBDMuseTwoStageTraining{ Name("SBDMuseTwoStageTraining"), Comment("If true, train two-stage model"), true };
        fhicl::Atom<std::string> SBDMallAtOnceModelFile{ Name("SBDMallAtOnceModelFile"), Comment("CSV filename for all-at-once model"), "" };
        fhicl::Atom<std::string> SBDMstage1ModelFile{ Name("SBDMstage1ModelFile"), Comment("CSV filename for stage-1 model"), "" };
        fhicl::Atom<std::string> SBDMstage2ModelFile{ Name("SBDMstage2ModelFile"), Comment("CSV filename for stage-2 model"), "" };
        fhicl::Atom<int> VirtualDetectorID{ Name("VirtualDetectorID"), Comment("Virtual detector ID to select"), 116 };
        fhicl::Atom<double> VDz0{ Name("VDz0"), Comment("z coordinate of the virtual detector"), 37700.39};
        fhicl::Atom<double> VDr{ Name("VDr"), Comment("VD radius"), 2000.0 };
        fhicl::Atom<int> pdgID{ Name("pdgID"), Comment("PDG ID to select"), 22 };
        fhicl::Atom<int> SBDMhidden{ Name("SBDMhidden"), Comment("Hidden layer size"), 128 };
        fhicl::Atom<int> SBDMlayers{ Name("SBDMlayers"), Comment("Number of layers"), 4 };
        fhicl::Atom<std::string> SBDMoptimizer{ Name("SBDMoptimizer"), Comment("Optimizer (SGD/ADAM)"), "ADAM" };
        fhicl::Atom<double> SBDMadamBeta1{ Name("SBDMadamBeta1"), Comment("Adam beta1"), 0.9 };
        fhicl::Atom<double> SBDMadamBeta2{ Name("SBDMadamBeta2"), Comment("Adam beta2"), 0.999 };
        fhicl::Atom<double> SBDMadamEps{ Name("SBDMadamEps"), Comment("Adam epsilon"), 1e-8 };
        fhicl::Atom<std::string> SBDMnoiseSchedule{ Name("SBDMnoiseSchedule"), Comment("Noise schedule (LINEAR/COSINE/LOGSIG)"), "COSINE" };
        fhicl::Atom<double> SBDMbetaMin{ Name("SBDMbetaMin"), Comment("Min beta (LINEAR)") , 1e-4 };
        fhicl::Atom<double> SBDMbetaMax{ Name("SBDMbetaMax"), Comment("Max beta (LINEAR)"), 0.02 };
        fhicl::Atom<double> SBDMcosineOffset{ Name("SBDMcosineOffset"), Comment("Cosine offset"), 0.008 };
        fhicl::Atom<double> SBDMlogSigMin{ Name("SBDMlogSigMin"), Comment("Minimum noise schedule parameter (for LOGSIG schedule)"), 1e-5 };
        fhicl::Atom<double> SBDMlogSigMax{ Name("SBDMlogSigMax"), Comment("Maximum noise schedule parameter (for LOGSIG schedule)"), 1.0 };
        fhicl::Atom<bool> SBDMepsPrediction{ Name("SBDMepsPrediction"), Comment("Whether predict eps (true) or the score (false)"), false };
        fhicl::Atom<double> SBDMlossWeightPower{ Name("SBDMlossWeightPower"), Comment("Power for weighting the loss function"), 2.0 };
        fhicl::Atom<int> SBDMbatchSize{ Name("SBDMbatchSize"), Comment("Batch size"), 32 };
        fhicl::Atom<double> SBDMgradientClip{ Name("SBDMgradientClip"), Comment("Gradient clip threshold"), 1.0 };
        fhicl::Atom<double> SBDMlearningRate{ Name("SBDMlearningRate"), Comment("Learning rate"), 1e-3 };
        fhicl::Atom<int> SBDMdiffusionSteps{ Name("SBDMdiffusionSteps"), Comment("Diffusion steps"), 200 };
        fhicl::Atom<int> SBDMtrainingSize{ Name("SBDMtrainingSize"), Comment("Size of the training data for the SBDM"), -1}; // -1 means use all available data
        fhicl::Atom<int> SBDMtrainingEpochs{ Name("SBDMtrainingEpochs"), Comment("Number of epochs"), 10 };
        fhicl::Sequence<int> SaveEpochs{ Name("SaveEpochs"), Comment("Epochs at which to save models"), std::vector<int>() };
        fhicl::Sequence<int> SBDMtrainingCurriculumEpochs{ Name("SBDMtrainingCurriculumEpochs"), Comment("Training curriculum epochs for each curriculum phase"), std::vector<int>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumLossWeightPower{ Name("SBDMtrainingCurriculumLossWeightPower"), Comment("Loss weight power for each curriculum phase"), std::vector<double>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumGradientClip{ Name("SBDMtrainingCurriculumGradientClip"), Comment("Gradient clip threshold for each curriculum phase"), std::vector<double>() };
        fhicl::Sequence<double> SBDMtrainingCurriculumLearningRate{ Name("SBDMtrainingCurriculumLearningRate"), Comment("Learning rate for each curriculum phase"), std::vector<double>() };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit VDResamplerTrainFromRoot(const Parameters& conf);
      void analyze(const art::Event&) override {}
      void endJob() override;
    private:
      art::RandomNumberGenerator::base_engine_t& engine;
      CLHEP::RandFlat randFlat_;
      CLHEP::RandGaussQ randGaussQ_;

      std::string inputRootFile;
      std::string treeName;
      bool useTwoStageTraining;
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
      std::vector<int> saveEpochs;

      int nPhase = 1;
      std::vector<int> trainingCurriculumEpochs = {};
      std::vector<double> trainingCurriculumLossWeightPower = {};
      std::vector<double> trainingCurriculumGradientClip = {};
      std::vector<double> trainingCurriculumLearningRate = {};
      std::vector<int> phaseBoundaries = {};

      // containers for normalization parameters
      double t_trans_mean = 0.0, t_trans_M2 = 0.0, t_trans_stdev = 0.0;
      double x_trans_mean = 0.0, x_trans_M2 = 0.0, x_trans_stdev = 0.0;
      double y_trans_mean = 0.0, y_trans_M2 = 0.0, y_trans_stdev = 0.0;
      double pr_t_mean = 0.0, pr_t_M2 = 0.0, pr_t_stdev = 0.0;
      double pphi_t_mean = 0.0, pphi_t_M2 = 0.0, pphi_t_stdev = 0.0;
      double pz_t_mean = 0.0, pz_t_M2 = 0.0, pz_t_stdev = 0.0;
      int nNorm = 0;
  };

  VDResamplerTrainFromRoot::VDResamplerTrainFromRoot(const Parameters& conf) :
    art::EDAnalyzer(conf),
    engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    randFlat_(engine),
    randGaussQ_(engine),
    inputRootFile(conf().InputRootFile()),
    treeName(conf().TreeName()),
    useTwoStageTraining(conf().SBDMuseTwoStageTraining()),
    SBDMallAtOnceModelFile(conf().SBDMallAtOnceModelFile()),
    SBDMstage1ModelFile(conf().SBDMstage1ModelFile()),
    SBDMstage2ModelFile(conf().SBDMstage2ModelFile()),
    VirtualDetectorID(conf().VirtualDetectorID()),
    VDz0(conf().VDz0()),
    VDr(conf().VDr()),
    pdgID(conf().pdgID()),
    trainingEpochs(conf().SBDMtrainingEpochs()),
    trainingSize(conf().SBDMtrainingSize()),
    saveEpochs(conf().SaveEpochs()),
    trainingCurriculumEpochs(conf().SBDMtrainingCurriculumEpochs()),
    trainingCurriculumLossWeightPower(conf().SBDMtrainingCurriculumLossWeightPower()),
    trainingCurriculumGradientClip(conf().SBDMtrainingCurriculumGradientClip()),
    trainingCurriculumLearningRate(conf().SBDMtrainingCurriculumLearningRate())
  {

    // Validate geometry configuration
    if (VDr <= 0.0) {
        throw cet::exception("VDResamplerTrainFromRoot")
            << "VDr must be positive (got " << VDr << "); "
            << "rho = r/VDr would produce inf/NaN in training data.";
    }
    if (!std::isfinite(VDz0)) {
        throw cet::exception("VDResamplerTrainFromRoot")
            << "VDz0 must be finite (got " << VDz0 << ").";
    }

    // validate curriculum parameters. If trainingCurriculumEpochs is not empty, the other vectors need
    // to be either the same length, or empty
    if (!trainingCurriculumEpochs.empty()) {
        nPhase = trainingCurriculumEpochs.size();
        if (nPhase == 1) {
            mf::LogWarning("VDResamplerTrainFromRoot")
                << "Only one curriculum phase. Curriculum training inputs will be ignored.";
        } else {
          std::stringstream ss;
          ss << "[Curriculum Training Schema] "<< std::endl;
          ss << nPhase << " phases."<< std::endl;
          trainingEpochs = std::accumulate(trainingCurriculumEpochs.begin(), trainingCurriculumEpochs.end(), 0);
          if (trainingCurriculumLossWeightPower.empty()) {
            for (int i = 0; i < nPhase; ++i) {
              trainingCurriculumLossWeightPower.push_back(conf().SBDMlossWeightPower());
            }
          } else if (static_cast<int>(trainingCurriculumLossWeightPower.size()) != nPhase) {
            throw cet::exception("VDResamplerTrainFromRoot")
                << "Inconsistent sizes for curriculum training parameters: SBDMtrainingCurriculumLossWeightPower.";
          }
          if (trainingCurriculumGradientClip.empty()) {
            for (int i = 0; i < nPhase; ++i) {
              trainingCurriculumGradientClip.push_back(conf().SBDMgradientClip());
            }
          } else if (static_cast<int>(trainingCurriculumGradientClip.size()) != nPhase) {
            throw cet::exception("VDResamplerTrainFromRoot")
                << "Inconsistent sizes for curriculum training parameters: SBDMtrainingCurriculumGradientClip.";
          }
          if (trainingCurriculumLearningRate.empty()) {
            for (int i = 0; i < nPhase; ++i) {
              trainingCurriculumLearningRate.push_back(conf().SBDMlearningRate());
            }
          } else if (static_cast<int>(trainingCurriculumLearningRate.size()) != nPhase) {
            throw cet::exception("VDResamplerTrainFromRoot")
                << "Inconsistent sizes for curriculum training parameters: SBDMtrainingCurriculumLearningRate.";
          }
          // Cache phase boundaries (ending epoch of each phase)
          phaseBoundaries.clear();
          int epochSum = 0;
          for (int i = 0; i < nPhase; ++i) {
            epochSum += trainingCurriculumEpochs[i];
            phaseBoundaries.push_back(epochSum);
          }
          ss << std::setw(10) << "Epochs"
             << std::setw(20) << "Loss Weight Power"
             << std::setw(20) << "Gradient Clip"
             << std::setw(20) << "Learning Rate"
             << std::endl;
          for (int i = 0; i < nPhase; ++i) {
            ss << std::setw(10) << trainingCurriculumEpochs[i]
               << std::setw(20) << trainingCurriculumLossWeightPower[i]
               << std::setw(20) << trainingCurriculumGradientClip[i]
               << std::setw(20) << trainingCurriculumLearningRate[i]
               << std::endl;
          }
          ss << "[End of Curriculum Training Schema]";
          mf::LogInfo("VDResamplerTrainFromRoot") << ss.str();
        }
    }

    // optimizer selection
    ScoreBasedDiffusionModel::OptimizerType opt;
    if (conf().SBDMoptimizer() == "SGD") {
        opt = ScoreBasedDiffusionModel::OptimizerType::SGD;
    } else {
        if (conf().SBDMoptimizer() != "ADAM") {
            mf::LogWarning("VDResamplerTrainFromRoot")
                << "Unrecognized SBDMoptimizer value \"" << conf().SBDMoptimizer()
                << "\"; falling back to ADAM.";
        }
        opt = ScoreBasedDiffusionModel::OptimizerType::ADAM;
    }

    // noise schedule
    ScoreBasedDiffusionModel::NoiseScheduleType sched;
    if (conf().SBDMnoiseSchedule() == "LINEAR") {
        sched = ScoreBasedDiffusionModel::NoiseScheduleType::LINEAR;
    } else if (conf().SBDMnoiseSchedule() == "LOGSIG") {
        sched = ScoreBasedDiffusionModel::NoiseScheduleType::LOGSIG;
    } else {
        if (conf().SBDMnoiseSchedule() != "COSINE") {
            mf::LogWarning("VDResamplerTrainFromRoot")
                << "Unrecognized SBDMnoiseSchedule value \"" << conf().SBDMnoiseSchedule()
                << "\"; falling back to COSINE.";
        }
        sched = ScoreBasedDiffusionModel::NoiseScheduleType::COSINE;
    }

    if (useTwoStageTraining) {
      stage1Model = std::make_unique<ScoreBasedDiffusionModel>(
      randFlat_, randGaussQ_, 3, 0, conf().SBDMhidden(), conf().SBDMlayers(), opt,
      conf().SBDMadamBeta1(), conf().SBDMadamBeta2(), conf().SBDMadamEps(), sched,
      conf().SBDMbetaMin(), conf().SBDMbetaMax(), conf().SBDMcosineOffset(),
      conf().SBDMlogSigMin(), conf().SBDMlogSigMax(), conf().SBDMepsPrediction(),
      nPhase > 1 ? trainingCurriculumLossWeightPower[0] : conf().SBDMlossWeightPower(),
      conf().SBDMbatchSize(),
      nPhase > 1 ? trainingCurriculumGradientClip[0] : conf().SBDMgradientClip(),
      nPhase > 1 ? trainingCurriculumLearningRate[0] : conf().SBDMlearningRate(),
      conf().SBDMdiffusionSteps()
      );
      stage2Model = std::make_unique<ScoreBasedDiffusionModel>(
      randFlat_, randGaussQ_, 3, 3, conf().SBDMhidden(), conf().SBDMlayers(), opt,
      conf().SBDMadamBeta1(), conf().SBDMadamBeta2(), conf().SBDMadamEps(), sched,
      conf().SBDMbetaMin(), conf().SBDMbetaMax(), conf().SBDMcosineOffset(),
      conf().SBDMlogSigMin(), conf().SBDMlogSigMax(), conf().SBDMepsPrediction(),
      nPhase > 1 ? trainingCurriculumLossWeightPower[0] : conf().SBDMlossWeightPower(),
      conf().SBDMbatchSize(),
      nPhase > 1 ? trainingCurriculumGradientClip[0] : conf().SBDMgradientClip(),
      nPhase > 1 ? trainingCurriculumLearningRate[0] : conf().SBDMlearningRate(),
      conf().SBDMdiffusionSteps()
      );
      if (trainingSize > 0) {
        stage1TrainingData.reserve(trainingSize);
        stage2TrainingData.reserve(trainingSize);
      } else {
        stage1TrainingData.reserve(1000);
        stage2TrainingData.reserve(1000);
      }
    } else {
      allAtOnceModel = std::make_unique<ScoreBasedDiffusionModel>(
      randFlat_, randGaussQ_, 6, 0, conf().SBDMhidden(), conf().SBDMlayers(), opt,
      conf().SBDMadamBeta1(), conf().SBDMadamBeta2(), conf().SBDMadamEps(), sched,
      conf().SBDMbetaMin(), conf().SBDMbetaMax(), conf().SBDMcosineOffset(),
      conf().SBDMlogSigMin(), conf().SBDMlogSigMax(), conf().SBDMepsPrediction(),
      nPhase > 1 ? trainingCurriculumLossWeightPower[0] : conf().SBDMlossWeightPower(),
      conf().SBDMbatchSize(),
      nPhase > 1 ? trainingCurriculumGradientClip[0] : conf().SBDMgradientClip(),
      nPhase > 1 ? trainingCurriculumLearningRate[0] : conf().SBDMlearningRate(),
      conf().SBDMdiffusionSteps()
      );
      if (trainingSize > 0) {
        allAtOnceTrainingData.reserve(trainingSize);
      } else {
        allAtOnceTrainingData.reserve(1000);
      }
    }
  }

  void VDResamplerTrainFromRoot::endJob() {
    // Open ROOT file and TTree
    TFile fin(inputRootFile.c_str(), "READ");
    if (!fin.IsOpen()) throw cet::exception("VDResamplerTrainFromRoot") << "Cannot open ROOT file: " << inputRootFile;
    TTree* ttree = dynamic_cast<TTree*>(fin.Get(treeName.c_str()));
    if (!ttree) throw cet::exception("VDResamplerTrainFromRoot") << "Cannot find TTree: " << treeName;

    // Set up branches
    double time, x, y, z, px, py, pz;
    int stepPdgId;
    ULong64_t virtualdetectorId;
    ttree->SetBranchAddress("time", &time);
    ttree->SetBranchAddress("x", &x);
    ttree->SetBranchAddress("y", &y);
    ttree->SetBranchAddress("z", &z);
    ttree->SetBranchAddress("px", &px);
    ttree->SetBranchAddress("py", &py);
    ttree->SetBranchAddress("pz", &pz);
    ttree->SetBranchAddress("pdgId", &stepPdgId);
    ttree->SetBranchAddress("virtualdetectorId", &virtualdetectorId);

    // Prepare training data
    for (Long64_t i = 0; i < ttree->GetEntries(); ++i) {
      ttree->GetEntry(i);

      if (virtualdetectorId != VirtualDetectorID || (stepPdgId != pdgID && pdgID != 0) || pz <= 0)
        continue; // Filter hits based on the virtual detector ID, particle type, and pz

      double x_trans, y_trans, t_trans, pr_t, pphi_t, pz_t;
      VDResampler::forwardTransformSample(x, y, z, time, px, py, pz,
                                          VDResampler::kX0, VDResampler::kY0, VDResampler::kT0, VDResampler::kTScale, VDResampler::kP0,
                                          VDr, VDz0,
                                          x_trans, y_trans, t_trans, pr_t, pphi_t, pz_t);
      if (nNorm < trainingSize || trainingSize <= 0) { // Wellford update
        nNorm++;
        double t_trans_delta = t_trans - t_trans_mean;
        t_trans_mean += t_trans_delta / static_cast<double>(nNorm);
        double t_trans_delta2 = t_trans - t_trans_mean;
        t_trans_M2 += t_trans_delta * t_trans_delta2;
        double x_trans_delta = x_trans - x_trans_mean;
        x_trans_mean += x_trans_delta / static_cast<double>(nNorm);
        double x_trans_delta2 = x_trans - x_trans_mean;
        x_trans_M2 += x_trans_delta * x_trans_delta2;
        double y_trans_delta = y_trans - y_trans_mean;
        y_trans_mean += y_trans_delta / static_cast<double>(nNorm);
        double y_trans_delta2 = y_trans - y_trans_mean;
        y_trans_M2 += y_trans_delta * y_trans_delta2;
        double pr_t_delta = pr_t - pr_t_mean;
        pr_t_mean += pr_t_delta / static_cast<double>(nNorm);
        double pr_t_delta2 = pr_t - pr_t_mean;
        pr_t_M2 += pr_t_delta * pr_t_delta2;
        double pphi_t_delta = pphi_t - pphi_t_mean;
        pphi_t_mean += pphi_t_delta / static_cast<double>(nNorm);
        double pphi_t_delta2 = pphi_t - pphi_t_mean;
        pphi_t_M2 += pphi_t_delta * pphi_t_delta2;
        double pz_t_delta = pz_t - pz_t_mean;
        pz_t_mean += pz_t_delta / static_cast<double>(nNorm);
        double pz_t_delta2 = pz_t - pz_t_mean;
        pz_t_M2 += pz_t_delta * pz_t_delta2;
      }
      if (useTwoStageTraining) {
        DiffusionTrainingSample s1; 
        s1.x = {t_trans, x_trans, y_trans};
        stage1TrainingData.push_back(std::move(s1));
        DiffusionTrainingSample s2; 
        s2.x = {pr_t, pphi_t, pz_t}; s2.cond = {t_trans, x_trans, y_trans};
        stage2TrainingData.push_back(std::move(s2));
      } else {
        DiffusionTrainingSample s; 
        s.x = {t_trans, x_trans, y_trans, pr_t, pphi_t, pz_t};
        allAtOnceTrainingData.push_back(std::move(s));
      }
    }
    fin.Close();

    t_trans_stdev = std::sqrt(t_trans_M2 / static_cast<double>(nNorm));
    x_trans_stdev = std::sqrt(x_trans_M2 / static_cast<double>(nNorm));
    y_trans_stdev = std::sqrt(y_trans_M2 / static_cast<double>(nNorm));
    pr_t_stdev = std::sqrt(pr_t_M2 / static_cast<double>(nNorm));
    pphi_t_stdev = std::sqrt(pphi_t_M2 / static_cast<double>(nNorm));
    pz_t_stdev = std::sqrt(pz_t_M2 / static_cast<double>(nNorm));

    if (useTwoStageTraining && (stage1TrainingData.empty() || stage2TrainingData.empty())) {
      mf::LogWarning("VDResamplerTrainFromRoot") << "No training data collected.";
      return;
    }
    if (!useTwoStageTraining && allAtOnceTrainingData.empty()) {
      mf::LogWarning("VDResamplerTrainFromRoot") << "No training data collected.";
      return;
    }

    // Training and saving at specified epochs
    if (useTwoStageTraining) {
      if(trainingSize > 0 && (int)stage1TrainingData.size() > trainingSize)
        stage1TrainingData.resize(trainingSize);
      if(trainingSize > 0 && (int)stage2TrainingData.size() > trainingSize)
        stage2TrainingData.resize(trainingSize);

      if (SBDMstage1ModelFile.empty() || SBDMstage2ModelFile.empty()) {
        throw cet::exception("VDResamplerTrainFromRoot") << "Two-stage training requires both SBDMstage1ModelFile and SBDMstage2ModelFile.";
      }

      // insert normalization information and normalize the training data
      std::vector<double> stage1Mean = {t_trans_mean, x_trans_mean, y_trans_mean};
      std::vector<double> stage1Stdev = {t_trans_stdev, x_trans_stdev, y_trans_stdev};
      stage1Model->normalizeData(stage1Mean, stage1Stdev, stage1TrainingData);
      std::vector<double> stage2Mean = {pr_t_mean, pphi_t_mean, pz_t_mean, t_trans_mean, x_trans_mean, y_trans_mean};
      std::vector<double> stage2Stdev = {pr_t_stdev, pphi_t_stdev, pz_t_stdev, t_trans_stdev, x_trans_stdev, y_trans_stdev};
      stage2Model->normalizeData(stage2Mean, stage2Stdev, stage2TrainingData);

      mf::LogInfo("VDResamplerTrainFromRoot")
          << "Training stage-1 diffusion model with " << stage1TrainingData.size()
          << " samples and " << trainingEpochs << " epochs...";
      std::string base = SBDMstage1ModelFile;
      if (base.size() > 4 && base.substr(base.size()-4) == ".csv") base = base.substr(0, base.size()-4);
      for (int e = 1; e <= trainingEpochs; ++e) {
        // Update training parameters at phase boundaries if curriculum is used
        if (nPhase > 1) {
          for (int phase = 0; phase < nPhase - 1; ++phase) {
            if (e == phaseBoundaries[phase] + 1) {
              double newLossWeightPower = trainingCurriculumLossWeightPower[phase + 1];
              double newGradientClip = trainingCurriculumGradientClip[phase + 1];
              double newLearningRate = trainingCurriculumLearningRate[phase + 1];
              mf::LogInfo("VDResamplerTrainFromRoot")
                << "Switching to phase " << (phase + 2)
                << ": lossWeightPower=" << newLossWeightPower
                << ", gradientClipThreshold=" << newGradientClip
                << ", learningRate=" << newLearningRate;
              stage1Model->updateLossWeightPower(newLossWeightPower);
              stage1Model->updateGradientClipThreshold(newGradientClip);
              stage1Model->updateLearningRate(newLearningRate);
            }
          }
        }
        stage1Model->train(stage1TrainingData, 1);
        if (std::find(saveEpochs.begin(), saveEpochs.end(), e) != saveEpochs.end()) {
          stage1Model->saveModel(base + ".epoch" + std::to_string(e) + ".csv");
        }
      }
      stage1Model->saveModel(SBDMstage1ModelFile); // Always save final model
      mf::LogInfo("VDResamplerTrainFromRoot") << "Stage-1 training completed";

      mf::LogInfo("VDResamplerTrainFromRoot")
          << "Training stage-2 diffusion model with " << stage2TrainingData.size()
          << " samples and " << trainingEpochs << " epochs...";
      base = SBDMstage2ModelFile;
      if (base.size() > 4 && base.substr(base.size()-4) == ".csv") base = base.substr(0, base.size()-4);
      for (int e = 1; e <= trainingEpochs; ++e) {
        // Update training parameters at phase boundaries if curriculum is used
        if (nPhase > 1) {
          for (int phase = 0; phase < nPhase - 1; ++phase) {
            if (e == phaseBoundaries[phase] + 1) {
              double newLossWeightPower = trainingCurriculumLossWeightPower[phase + 1];
              double newGradientClip = trainingCurriculumGradientClip[phase + 1];
              double newLearningRate = trainingCurriculumLearningRate[phase + 1];
              mf::LogInfo("VDResamplerTrainFromRoot")
                << "Switching to phase " << (phase + 2)
                << ": lossWeightPower=" << newLossWeightPower
                << ", gradientClipThreshold=" << newGradientClip
                << ", learningRate=" << newLearningRate;
              stage2Model->updateLossWeightPower(newLossWeightPower);
              stage2Model->updateGradientClipThreshold(newGradientClip);
              stage2Model->updateLearningRate(newLearningRate);
            }
          }
        }
        stage2Model->train(stage2TrainingData, 1);
        if (std::find(saveEpochs.begin(), saveEpochs.end(), e) != saveEpochs.end()) {
          stage2Model->saveModel(base + ".epoch" + std::to_string(e) + ".csv");
        }
      }
      stage2Model->saveModel(SBDMstage2ModelFile);
      mf::LogInfo("VDResamplerTrainFromRoot") << "Stage-2 training completed";

    } else {
      if(trainingSize > 0 && (int)allAtOnceTrainingData.size() > trainingSize)
        allAtOnceTrainingData.resize(trainingSize);

      if (SBDMallAtOnceModelFile.empty()) {
        throw cet::exception("VDResamplerTrainFromRoot") << "All-at-once training requires SBDMallAtOnceModelFile.";
      }

      // insert normalization information and normalize the training data
      std::vector<double> allAtOnceMean = {t_trans_mean, x_trans_mean, y_trans_mean, pr_t_mean, pphi_t_mean, pz_t_mean};
      std::vector<double> allAtOnceStdev = {t_trans_stdev, x_trans_stdev, y_trans_stdev, pr_t_stdev, pphi_t_stdev, pz_t_stdev};
      allAtOnceModel->normalizeData(allAtOnceMean, allAtOnceStdev, allAtOnceTrainingData);

      mf::LogInfo("VDResamplerTrainFromRoot")
          << "Training all-at-once diffusion model with " << allAtOnceTrainingData.size()
          << " samples and " << trainingEpochs << " epochs...";
      std::string base = SBDMallAtOnceModelFile;
      if (base.size() > 4 && base.substr(base.size()-4) == ".csv") base = base.substr(0, base.size()-4);
      for (int e = 1; e <= trainingEpochs; ++e) {
        // Update training parameters at phase boundaries if curriculum is used
        if (nPhase > 1) {
          for (int phase = 0; phase < nPhase - 1; ++phase) {
            if (e == phaseBoundaries[phase] + 1) {
              double newLossWeightPower = trainingCurriculumLossWeightPower[phase + 1];
              double newGradientClip = trainingCurriculumGradientClip[phase + 1];
              double newLearningRate = trainingCurriculumLearningRate[phase + 1];
              mf::LogInfo("VDResamplerTrainFromRoot")
                << "Switching to phase " << (phase + 2)
                << ": lossWeightPower=" << newLossWeightPower
                << ", gradientClipThreshold=" << newGradientClip
                << ", learningRate=" << newLearningRate;
              allAtOnceModel->updateLossWeightPower(newLossWeightPower);
              allAtOnceModel->updateGradientClipThreshold(newGradientClip);
              allAtOnceModel->updateLearningRate(newLearningRate);
            }
          }
        }
        allAtOnceModel->train(allAtOnceTrainingData, 1);
        if (std::find(saveEpochs.begin(), saveEpochs.end(), e) != saveEpochs.end()) {
          allAtOnceModel->saveModel(base + ".epoch" + std::to_string(e) + ".csv");
        }
      }
      allAtOnceModel->saveModel(SBDMallAtOnceModelFile);
      mf::LogInfo("VDResamplerTrainFromRoot") << "All-at-once training completed";
    }
    return;
  }
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::VDResamplerTrainFromRoot)
