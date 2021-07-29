// Generates electrons in a flat shaped spectrum, size determined by FCL params these will be attached to a mu- in
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
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataTable.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"

#include <TTree.h>
namespace mu2e {

  //================================================================
  class FlatMuonDaughterGenerator : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<double> startMom{Name("startMom"),0};
      fhicl::Atom<double> endMom{Name("endMom"),105};
      fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),Comment("A SimParticleCollection with input stopped muons.")};
      fhicl::Atom<std::string> stoppingTargetMaterial{Name("stoppingTargetMaterial"),Comment("material"),"Al" };
      fhicl::Atom<unsigned> verbosity{Name("verbosity"),0};
      fhicl::Atom<bool> makeHistograms{Name("makeHistograms"),false};
      fhicl::Atom<std::string> processcode{Name("processcode")};
      fhicl::Atom<int> pdgId{Name("pdgId")};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit FlatMuonDaughterGenerator(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:
    Float_t pmag_gen;
    Float_t time_gen;
    Float_t x_gen;
    Float_t y_gen;
    Float_t z_gen;
    TTree*  _Ntup;
    const PDGCode::type electronId_ = PDGCode::e_minus;
    double electronMass_;
    double startMom_;
    double endMom_;
    double muonLifeTime_;

    art::ProductToken<SimParticleCollection> const simsToken_;

    unsigned verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randFlat_;
    CLHEP::RandExponential randExp_;
    RandomUnitSphere   randomUnitSphere_;
    bool makeHistograms_;
    std::string processcode_;
    int pdgId_;
  };

  //================================================================
  FlatMuonDaughterGenerator::FlatMuonDaughterGenerator(const Parameters& conf)
    : EDProducer{conf}
    , electronMass_(GlobalConstantsHandle<ParticleDataTable>()->particle(electronId_).ref().mass().value())
    , startMom_(conf().startMom())
    , endMom_(conf().endMom())
    , muonLifeTime_{GlobalConstantsHandle<PhysicsParams>()->getDecayTime(conf().stoppingTargetMaterial())}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randFlat_{eng_}
    , randExp_{eng_}
    , randomUnitSphere_{eng_}
    , makeHistograms_(conf().makeHistograms())
    , processcode_(conf().processcode())
    , pdgId_(conf().pdgId())
  {
    produces<mu2e::StageParticleCollection>();


    if(makeHistograms_){
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "FlatGun" );
      _Ntup  = tfs->make<TTree>("GenTree", "GenTree");
	    _Ntup->Branch("pmag_gen", &pmag_gen , "pmag_gen/F");
	    _Ntup->Branch("time_gen", &time_gen, "time_gen/F");
	    _Ntup->Branch("x_gen", &x_gen, "x_gen/F");
	    _Ntup->Branch("y_gen", &y_gen, "y_gen/F");
	    _Ntup->Branch("z_gen", &z_gen, "z_gen/F");
	    }
  }

  //================================================================
  void FlatMuonDaughterGenerator::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto mus = stoppedMuMinusList(simh);

    if(mus.empty()) {
      throw   cet::exception("BADINPUT")
        <<"FlatMuonDaughterGenerator::produce(): no suitable stopped mu- in the input SimParticleCollection\n";

    }
    
    const auto mustop = mus.at(eng_.operator unsigned int() % mus.size());
    double randomMom = randFlat_.fire(startMom_, endMom_);
    double randomE = sqrt(electronMass_*electronMass_ + randomMom*randomMom);
    double time = mustop->endGlobalTime() + randExp_.fire(muonLifeTime_);
    ProcessCode pcode_ = ProcessCode::findByName(processcode_);
    PDGCode::type pid = static_cast<PDGCode::type>(pdgId_);
    output->emplace_back(mustop,
                         pcode_,
                         pid,
                         mustop->endPosition(),
                         CLHEP::HepLorentzVector{randomUnitSphere_.fire(randomMom), randomE},
                         time
                         );

    if(makeHistograms_){
      pmag_gen = randomMom;
      time_gen = time;
      x_gen = mustop->endPosition().x();
      y_gen = mustop->endPosition().y();
      z_gen = mustop->endPosition().z();
      _Ntup->Fill();
    }
    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FlatMuonDaughterGenerator);
