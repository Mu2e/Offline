// VDResamplerGenerateFromModel_module.cc
// This module takes either a single all-at-once model file or a pair of stage-1/stage-2
// model files and produces new virtual-detector-like GenParticles from the trained
// score-based diffusion model(s).
// Yongyi Wu, Mar. 2026

// stdlib includes
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

#include "Offline/MachineLearningTools/inc/ScoreBasedDiffusionModel.hh"

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/SeedService/inc/SeedService.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

namespace mu2e {

  namespace {
    // Utility function to extract PDG ID from the model filename, which is expected to contain a substring like "pdg[optional m][digits]", e.g. "pdg13" for muons, "pdgm13" for muon pluses.
    int loadPDGIdFromFileName(const std::string& fileName) {
      const size_t pdgPos = fileName.find("pdg");
      if (pdgPos == std::string::npos) {
        throw cet::exception("VDResamplerGenerateFromModel")
          << "Cannot infer PDG ID from model filename: " << fileName;
      }
      size_t pos = pdgPos + 3;
      bool negative = false;
      if (pos < fileName.size() && (fileName[pos] == 'm' || fileName[pos] == '-')) {
        negative = true;
        ++pos;
      }
      const size_t startDigits = pos;
      while (pos < fileName.size() && std::isdigit(static_cast<unsigned char>(fileName[pos]))) {
        ++pos;
      }
      if (startDigits == pos) {
        throw cet::exception("VDResamplerGenerateFromModel")
          << "Cannot infer PDG ID from model filename: " << fileName;
      }
      const int magnitude = std::stoi(fileName.substr(startDigits, pos - startDigits));
      return negative ? -magnitude : magnitude;
    }
  }

  class VDResamplerGenerateFromModel : public art::EDProducer {
    public:
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      struct Config {
        fhicl::Atom<bool> useTwoStageModel{
          Name("useTwoStageModel"),
          Comment("If true, use stage-1/stage-2 models. If false, use the all-at-once 6D model."),
          true
        };
        fhicl::Atom<std::string> stage1ModelFile{
          Name("stage1ModelFile"),
          Comment("CSV filename for the stage-1 model parameters"),
          ""
        };
        fhicl::Atom<std::string> stage2ModelFile{
          Name("stage2ModelFile"),
          Comment("CSV filename for the stage-2 model parameters"),
          ""
        };
        fhicl::Atom<std::string> allAtOnceModelFile{
          Name("allAtOnceModelFile"),
          Comment("CSV filename for the all-at-once 6D model parameters"),
          ""
        };
        fhicl::Atom<bool> useHeun{
          Name("useHeun"),
          Comment("If true, use Heun's method for reverse diffusion. Otherwise use Euler."),
          true
        };
        fhicl::Atom<int> diffusionSteps{
          Name("diffusionSteps"),
          Comment("Number of reverse-diffusion steps used for sampling"),
          200
        };
        fhicl::Atom<double> VDz0{
          Name("VDz0"),
          Comment("Nominal z coordinate of the virtual detector"),
          37700.39
        };
        fhicl::Atom<double> VDr{
          Name("VDr"),
          Comment("Virtual detector radius used in the training transform"),
          2000.0
        };
        fhicl::Atom<bool> doROOTDump{
          Name("doROOTDump"),
          Comment("Whether to dump generated samples into a ROOT tree"),
          false
        };
      };

      using Parameters = art::EDProducer::Table<Config>;

      explicit VDResamplerGenerateFromModel(const Parameters& conf);
      void produce(art::Event& event) override;

    private:
      art::RandomNumberGenerator::base_engine_t& engine_;
      bool useTwoStageModel_;
      std::string stage1ModelFile_;
      std::string stage2ModelFile_;
      std::string allAtOnceModelFile_;
      bool useHeun_;
      int diffusionSteps_;
      double VDz0_;
      double VDr_;
      bool doROOTDump_;
      GlobalConstantsHandle<ParticleDataList> pdt_;

      std::unique_ptr<ScoreBasedDiffusionModel> allAtOnceModel_;
      std::unique_ptr<ScoreBasedDiffusionModel> stage1Model_;
      std::unique_ptr<ScoreBasedDiffusionModel> stage2Model_;

      int pdgId_ = 0;

      // Detector-center parameters used in the same transform as training.
      double x0_ = -3904.0;
      double y0_ = 0.0;
      double t0_ = 1700.0;   // ns
      double tScale_ = 1.0;
      double p0_ = 1.0;

      // Variables for optional ROOT dump.
      TTree* outTree_ = nullptr;
      double x_gen_ = 0.0;
      double y_gen_ = 0.0;
      double z_gen_ = 0.0;
      double t_gen_ = 0.0;
      double px_gen_ = 0.0;
      double py_gen_ = 0.0;
      double pz_gen_ = 0.0;
      double mass_gen_ = 0.0;
      double E_gen_ = 0.0;
    };

  VDResamplerGenerateFromModel::VDResamplerGenerateFromModel(const Parameters& conf)
    : art::EDProducer{conf},
      engine_(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
      useTwoStageModel_(conf().useTwoStageModel()),
      stage1ModelFile_(conf().stage1ModelFile()),
      stage2ModelFile_(conf().stage2ModelFile()),
      allAtOnceModelFile_(conf().allAtOnceModelFile()),
      useHeun_(conf().useHeun()),
      diffusionSteps_(conf().diffusionSteps()),
      VDz0_(conf().VDz0()),
      VDr_(conf().VDr()),
      doROOTDump_(conf().doROOTDump()) {

    produces<GenParticleCollection>();

    if (useTwoStageModel_) {
      if (stage1ModelFile_.empty() || stage2ModelFile_.empty()) {
        throw cet::exception("VDResamplerGenerateFromModel")
          << "Two-stage generation requires both stage1ModelFile and stage2ModelFile.";
      }

      stage1Model_ = std::make_unique<ScoreBasedDiffusionModel>(
        ScoreBasedDiffusionModel::loadModel(engine_, stage1ModelFile_)
      );
      stage2Model_ = std::make_unique<ScoreBasedDiffusionModel>(
        ScoreBasedDiffusionModel::loadModel(engine_, stage2ModelFile_)
      );
      pdgId_ = loadPDGIdFromFileName(stage1ModelFile_);
    } else {
      if (allAtOnceModelFile_.empty()) {
        throw cet::exception("VDResamplerGenerateFromModel")
          << "All-at-once generation requires allAtOnceModelFile.";
      }

      allAtOnceModel_ = std::make_unique<ScoreBasedDiffusionModel>(
        ScoreBasedDiffusionModel::loadModel(engine_, allAtOnceModelFile_)
      );
      pdgId_ = loadPDGIdFromFileName(allAtOnceModelFile_);
    }

    z_gen_ = VDz0_;

    if (doROOTDump_) {
      art::ServiceHandle<art::TFileService> tfs;
      outTree_ = tfs->make<TTree>("ttree", "Generated samples from VD resampler model");
      outTree_->Branch("pdgId", &pdgId_, "pdgId/I");
      outTree_->Branch("x", &x_gen_, "x/D");
      outTree_->Branch("y", &y_gen_, "y/D");
      outTree_->Branch("z", &z_gen_, "z/D");
      outTree_->Branch("time", &t_gen_, "time/D");
      outTree_->Branch("px", &px_gen_, "px/D");
      outTree_->Branch("py", &py_gen_, "py/D");
      outTree_->Branch("pz", &pz_gen_, "pz/D");
      outTree_->Branch("E", &E_gen_, "E/D");
    }
  }

  void VDResamplerGenerateFromModel::produce(art::Event& event) {
    auto output = std::make_unique<GenParticleCollection>();

    double x_trans = 0.0;
    double y_trans = 0.0;
    double t_trans = 0.0;
    double pr_t = 0.0;
    double pphi_t = 0.0;
    double pz_t = 0.0;

    // Generate a new sample using the loaded model(s)
    // note the values here are transformed and need to be inverted back to the original coordinates after sampling.
    if (useTwoStageModel_) {
      const std::vector<double> stage1Sample = stage1Model_->generateSample({}, useHeun_, diffusionSteps_);
      if (stage1Sample.size() != 3u) {
        throw cet::exception("VDResamplerGenerateFromModel")
          << "Stage-1 model returned " << stage1Sample.size() << " values, expected 3.";
      }

      t_trans = stage1Sample[0];
      x_trans = stage1Sample[1];
      y_trans = stage1Sample[2];

      const std::vector<double> stage2Condition = {t_trans, x_trans, y_trans};
      const std::vector<double> stage2Sample = stage2Model_->generateSample(stage2Condition, useHeun_, diffusionSteps_);
      if (stage2Sample.size() != 3u) {
        throw cet::exception("VDResamplerGenerateFromModel")
          << "Stage-2 model returned " << stage2Sample.size() << " values, expected 3.";
      }

      pr_t = stage2Sample[0];
      pphi_t = stage2Sample[1];
      pz_t = stage2Sample[2];
    } else {
      const std::vector<double> sample = allAtOnceModel_->generateSample({}, useHeun_, diffusionSteps_);
      if (sample.size() != 6u) {
        throw cet::exception("VDResamplerGenerateFromModel")
          << "All-at-once model returned " << sample.size() << " values, expected 6.";
      }

      t_trans = sample[0];
      x_trans = sample[1];
      y_trans = sample[2];
      pr_t = sample[3];
      pphi_t = sample[4];
      pz_t = sample[5];
    }

    // Recover detector coordinates from the transformed (x', y').
    const double u = std::sqrt(x_trans * x_trans + y_trans * y_trans);
    const double theta = std::atan2(y_trans, x_trans);
    const double rho = std::tanh(u);
    const double r = rho * VDr_;
    const double dx = r * std::cos(theta);
    const double dy = r * std::sin(theta);
    x_gen_ = dx + x0_;
    y_gen_ = dy + y0_;
    z_gen_ = VDz0_;

    // Invert the time transform.
    t_gen_ = t0_ * std::exp(t_trans * tScale_);

    // Invert the momentum transform.
    const double pr = p0_ * std::sinh(pr_t);
    const double pphi = p0_ * std::sinh(pphi_t);
    pz_gen_ = p0_ * std::sinh(pz_t);

    if (r > 1e-6) {
      const double rx = dx / r;
      const double ry = dy / r;
      const double phix = -ry;
      const double phiy = rx;
      px_gen_ = pr * rx + pphi * phix;
      py_gen_ = pr * ry + pphi * phiy;
    } else {
      px_gen_ = pr;
      py_gen_ = pphi;
    }

    mass_gen_ = pdt_->particle(pdgId_).mass();
    const CLHEP::Hep3Vector momParticle(px_gen_, py_gen_, pz_gen_);
    const CLHEP::Hep3Vector posParticle(x_gen_, y_gen_, z_gen_);
    const double eTotal = std::sqrt(momParticle.mag2() + mass_gen_ * mass_gen_);
    E_gen_ = eTotal - mass_gen_;
    const CLHEP::HepLorentzVector fourMomParticle(eTotal, momParticle);

    output->emplace_back(
      PDGCode::type(pdgId_),
      GenId::STMDownStreamGenTool,
      posParticle,
      fourMomParticle,
      t_gen_
    );

    event.put(std::move(output));

    if (doROOTDump_) {
      outTree_->Fill();
    }
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::VDResamplerGenerateFromModel)
