// This module will be used to make electrons or positrons from any process originating from mu- or mu+ stopped in any target material
//
// Sophie Middleton, 2021
// S. Huang updates for mu+ stops

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandExponential.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
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
                  "Only aluminum (Al) is supported, emisson spectra for other materials are not implemented.\n"),
          "Al" };
      //  fhicl::Atom<int> muoncharge{Name("muoncharge"), Comment("charge = -1: mu-; charge = +1: mu+.\n" ),-1};  
      fhicl::Atom<int> pdgId{Name("pdgId"),Comment("pdg id of daughter particle"),11};
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
    //   int muoncharge_;
    int pdgId_;


    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;
  
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
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
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

    
  }

  //================================================================
  void SingleProcessGenerator::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    // int eventid=event.id().event();  
    // int count_particle=0;
    auto mus_temp = stoppedMuMinusList(simh); //default, negative muon stops
 
    if (pdgId_==-11){
      mus_temp = stoppedMuPlusList(simh); //positive muon stops
      muonLifeTime_=0; //decay time already included for stopped muon(+)
    }
   
    const auto mus=mus_temp;

     for(const auto& mustop: mus) {
    //  count_particle++;
    
      const double time = mustop->endGlobalTime() + randExp_.fire(muonLifeTime_);

      addParticles(output.get(), mustop, time, Generator_.get());
      //  printf("the event number is %d\n",eventid);
      //   if(count_particle>1)printf("number of particle is %d\t in event %d\n",count_particle,eventid);
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

DEFINE_ART_MODULE(mu2e::SingleProcessGenerator);
