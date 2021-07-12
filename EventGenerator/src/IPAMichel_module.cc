// S Middleton 2021

#include <iostream>
#include <string>
#include <cmath>
#include <memory>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/StageParticle.hh"
#include "Mu2eUtilities/inc/simParticleList.hh"

namespace mu2e {

  //================================================================
  class IPAMichel : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),
          Comment("A SimParticleCollection with input stopped muons.")};
      fhicl::Atom<unsigned> verbosity{Name("verbosity"),0};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit IPAMichel(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:
    const PDGCode::type electronId_ = PDGCode::e_minus;
    double electronMass_;
    double endPointEnergy_;
    double endPointMomentum_;
    double muonLifeTime_;

    art::ProductToken<SimParticleCollection> const simsToken_;

    unsigned verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;
    RandomUnitSphere   randomUnitSphere_;
  };

  //================================================================
  IPAMichel::IPAMichel(const Parameters& conf)
    : EDProducer{conf}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
    , randomUnitSphere_{eng_}
  {
    produces<mu2e::StageParticleCollection>();
    
  }

  //================================================================
  void IPAMichel::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto mus = stoppedMuMinusList(simh);

    if(mus.empty()) {
      throw   cet::exception("BADINPUT")
        <<"IPAMichel::produce(): no suitable stopped mu- in the input SimParticleCollection\n";

    }

    const auto mustop = mus.at(eng_.operator unsigned int() % mus.size());

    output->emplace_back(mustop,
                         ProcessCode::mu2eMuonDecayAtRest, //process code for DIO (TODO - unique code for IPAMichel?)
                         PDGCode::e_minus,
                         mustop->endPosition(), //final position of stop
                         CLHEP::HepLorentzVector{randomUnitSphere_.fire(endPointMomentum_), endPointEnergy_}, //random choice of momentum/energy 
                         mustop->endGlobalTime() + randExp_.fire(muonLifeTime_) //random time choice
                         );


    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::IPAMichel);
