// Generates electron at the endpoint energy that will be attached to a mu- in
// the input SimParticleCollection.
// This module throws an exception if no suitable muon is found.
//
// Andrei Gaponenko, 2021

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

namespace mu2e {

  //================================================================
  class CeEndpoint : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),
          Comment("A SimParticleCollection with input stopped muons.")};
      fhicl::Atom<std::string> stoppingTargetMaterial{
        Name("stoppingTargetMaterial"),
          Comment("Material determines endpoint energy and muon life time.  Material must be known to the GlobalConstantsService."),
          "Al" };
      fhicl::Atom<unsigned> verbosity{Name("verbosity"),0};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit CeEndpoint(const Parameters& conf);

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
  CeEndpoint::CeEndpoint(const Parameters& conf)
    : EDProducer{conf}
    , electronMass_(GlobalConstantsHandle<ParticleDataTable>()->particle(electronId_).ref().mass().value())
    , endPointEnergy_{GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy(conf().stoppingTargetMaterial())}
    , endPointMomentum_{ endPointEnergy_*sqrt(1 - std::pow(electronMass_/endPointEnergy_,2)) }
    , muonLifeTime_{GlobalConstantsHandle<PhysicsParams>()->getDecayTime(conf().stoppingTargetMaterial())}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
    , randomUnitSphere_{eng_}
  {
    produces<mu2e::StageParticleCollection>();
    if(verbosity_ > 0) {
      mf::LogInfo log("CeEndpoint");
      log<<"stoppingTargetMaterial = "<<conf().stoppingTargetMaterial()
         <<", endpoint energy = "<<endPointEnergy_
         <<", muon lifetime = "<<muonLifeTime_
         <<std::endl;
    }
  }

  //================================================================
  void CeEndpoint::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);

    const auto& sims = *simh;
    if(!sims.empty()) {
      const auto stage = sims.back().second.simStage();

      // We want to find a single stopped muon.   Even if an event contains multiple
      // mu stops, we should not create multiple conversion electrons per event.
      art::Ptr<SimParticle> mustop;

      for(auto i = sims.rbegin(); i != sims.rend(); ++i) {
        const auto& inpart = i->second;
        if((inpart.simStage() == stage) &&
           (inpart.pdgId() == PDGCode::mu_minus) &&
           (inpart.stoppingCode() == ProcessCode::muMinusCaptureAtRest)) // G4 sets this code for both decay and capture cases
          {
            mustop = art::Ptr<SimParticle>(simh, i->first.asInt());
            break;
          }
      }

      if(mustop) {
        output->emplace_back(mustop,
                             ProcessCode::muMinusConversionAtRest,
                             PDGCode::e_minus,
                             mustop->endPosition(),
                             CLHEP::HepLorentzVector{randomUnitSphere_.fire(endPointMomentum_), endPointEnergy_},
                             mustop->endGlobalTime() + randExp_.fire(muonLifeTime_)
                             );
      }
    }

    const bool passed = !output->empty();
    event.put(std::move(output));
    if(!passed) {
      throw   cet::exception("BADINPUT")
        <<"CeEndpoint::produce(): no suitable stopped mu- in the input SimParticleCollection\n";

    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CeEndpoint);
