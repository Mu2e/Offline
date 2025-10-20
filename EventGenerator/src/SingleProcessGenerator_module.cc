// This module will be used to make electrons or positrons from any process originating from mu- stopped in any target material
//
// Sophie Middleton, 2021


#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"

#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

namespace mu2e {
  //================================================================
  class SingleProcessGenerator : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),
          Comment("A SimParticleCollection with muons stopping in the target.")};
      fhicl::Atom<std::string> stoppingTargetMaterial{
        Name("stoppingTargetMaterial"),
          Comment("Material determines muon life time, capture fraction, and particle spectra.\n"
                  "Only aluminum (Al) is supported, emisson spectra for other materials are not implemented.\n"),"Al" };
      fhicl::Atom<int> pdgId{Name("pdgId"),Comment("pdg id of daughter particle"),PDGCode::e_minus};
      fhicl::Atom<bool> selectOne{Name("selectOne"),Comment("Select one input muon from the stopped muon list"), false};

      fhicl::DelegatedParameter decayProducts{Name("decayProducts"), Comment("spectrum (and variables) to be generated")};
      fhicl::Atom<unsigned> verbosity{Name("verbosity"),0};

    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit SingleProcessGenerator(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:
    double muonLifeTime_;
    art::ProductToken<SimParticleCollection> const simsToken_;
    unsigned verbosity_;
    int pdgId_;
    bool selectOne_;


    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;
    CLHEP::RandFlat randFlat_;

    std::unique_ptr<ParticleGeneratorTool> Generator_;

    void addParticles(StageParticleCollection* output, art::Ptr<SimParticle> mustop, double time, ParticleGeneratorTool* gen);
  };

  //================================================================
  SingleProcessGenerator::SingleProcessGenerator(const Parameters& conf)
    : EDProducer{conf}
    , muonLifeTime_{GlobalConstantsHandle<PhysicsParams>()->getDecayTime(conf().stoppingTargetMaterial())}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , pdgId_{conf().pdgId()}
    , selectOne_{conf().selectOne()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
    , randFlat_{eng_}
  {
    produces<mu2e::StageParticleCollection>();
    if(verbosity_ > 0) {
      mf::LogInfo log("SingleProcessGenerator");
      log<<"stoppingTargetMaterial = "<<conf().stoppingTargetMaterial()
         <<", muon lifetime = "<<muonLifeTime_
         <<std::endl;
    }

    const auto pset = conf().decayProducts.get<fhicl::ParameterSet>();
    Generator_ = art::make_tool<ParticleGeneratorTool>(pset);
    Generator_->finishInitialization(eng_, conf().stoppingTargetMaterial());

    if(pdgId_==PDGCode::e_plus) {
      muonLifeTime_=0; //decay time already included for stopped muon(+) FIXME!!!
    }
  }

  //================================================================
  void SingleProcessGenerator::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto mus  = (pdgId_==PDGCode::e_minus) ? stoppedMuMinusList(simh) : stoppedMuPlusList(simh);

    const bool   select_one  = selectOne_ && !mus.empty();
    const size_t first_index = (select_one) ? randFlat_.fireInt(mus.size()) : 0;
    const size_t last_index  = (select_one) ? first_index + 1               : mus.size(); // last index + 1

    for(size_t index = first_index; index < last_index; ++index) {
      const auto& mustop = mus[index];
      const double time  = mustop->endGlobalTime() + randExp_.fire(muonLifeTime_);
      addParticles(output.get(), mustop, time, Generator_.get());
    }

    event.put(std::move(output));
  }

  //================================================================
  void SingleProcessGenerator::addParticles(StageParticleCollection* output,
                            art::Ptr<SimParticle> mustop,
                            double time,
                            ParticleGeneratorTool* gen)
  {
    auto daughters = gen->generate();
    for(const auto& d: daughters) {

      output->emplace_back(mustop,
                           d.creationCode,
                           d.pdgId,
                           mustop->endPosition(),
                           d.fourmom,
                           time
                           );

    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SingleProcessGenerator)
