// Generates particles in a flat shaped spectrum, size determined by FCL params these will be attached to a mu- in
// the input SimParticleCollection.
// This module throws an exception if no suitable muon is found.
//
// S Middleton, 2021

#include <iostream>
#include <string>
#include <cmath>
#include <memory>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"

#include <TTree.h>
namespace mu2e {

  //================================================================
  class FlatGeneratorFromStoppedMuon : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<double> startMom{Name("startMom"),0};
      fhicl::Atom<double> endMom{Name("endMom"),105};
      fhicl::Atom<double> czMin{Name("czMin"), Comment("Minimum cos(theta_z) to generate in"), -1.};
      fhicl::Atom<double> czMax{Name("czMax"), Comment("Maximum cos(theta_z) to generate in"),  1.};
      fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),Comment("A SimParticleCollection with input stopped muons.")};
      fhicl::Atom<std::string> stoppingTargetMaterial{Name("stoppingTargetMaterial"),Comment("material")};
      fhicl::Atom<unsigned> verbosity{Name("verbosity")};
      fhicl::Atom<int> pdgId{Name("pdgId"),Comment("pdg id of daughter particle")};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit FlatGeneratorFromStoppedMuon(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:

    double particleMass_;
    double startMom_;
    double endMom_;
    double muonLifeTime_;

    art::ProductToken<SimParticleCollection> const simsToken_;

    unsigned verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randFlat_;
    CLHEP::RandExponential randExp_;
    RandomUnitSphere   randomUnitSphere_;
    ProcessCode process;
    int pdgId_;
    PDGCode::type pid;
  };

  //================================================================
  FlatGeneratorFromStoppedMuon::FlatGeneratorFromStoppedMuon(const Parameters& conf)
    : EDProducer{conf}
    , particleMass_(GlobalConstantsHandle<ParticleDataList>()->particle(static_cast<PDGCode::type>(conf().pdgId())).mass())
    , startMom_(conf().startMom())
    , endMom_(conf().endMom())
    , muonLifeTime_{GlobalConstantsHandle<PhysicsParams>()->getDecayTime(conf().stoppingTargetMaterial())}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randFlat_{eng_}
    , randExp_{eng_}
    , randomUnitSphere_{eng_, conf().czMin(), conf().czMax()}
    , pdgId_(conf().pdgId())

  {
    produces<mu2e::StageParticleCollection>();
    pid = static_cast<PDGCode::type>(pdgId_);

    if (pid == PDGCode::e_minus) { process = ProcessCode::mu2eFlateMinus; }
    else if (pid == PDGCode::e_plus) { process = ProcessCode::mu2eFlatePlus; }
    else if (pid == PDGCode::gamma) { process = ProcessCode::mu2eFlatPhoton; }
    else if (pid == PDGCode::mu_minus) { process = ProcessCode::mu2eFlatMuMinus; }
    else {
      throw   cet::exception("BADINPUT")
        <<"FlatGeneratorFromStoppedMuon::produce(): No process associated with chosen PDG id\n";
    }

  }

  //================================================================
  void FlatGeneratorFromStoppedMuon::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto mus = stoppedMuMinusList(simh);

    if(mus.empty()) {
      throw   cet::exception("BADINPUT")
        <<"FlatGeneratorFromStoppedMuon::produce(): no suitable stopped mu- in the input SimParticleCollection\n";

    }

    const auto mustop = mus.at(eng_.operator unsigned int() % mus.size());
    double randomMom = randFlat_.fire(startMom_, endMom_);
    double randomE = sqrt(particleMass_*particleMass_ + randomMom*randomMom);
    double time = mustop->endGlobalTime() + randExp_.fire(muonLifeTime_);


    output->emplace_back(mustop,
                         process,
                         pid,
                         mustop->endPosition(),
                         CLHEP::HepLorentzVector{randomUnitSphere_.fire(randomMom), randomE},
                         time
                         );

    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FlatGeneratorFromStoppedMuon)
