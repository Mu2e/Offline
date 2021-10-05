// Generates electron with Leading Log spectrum that will be attached to a mu- in
// the input SimParticleCollection.
// This module throws an exception if no suitable muon is found.
//
// Sophie Middleton, 2021

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
#include "fhiclcpp/types/DelegatedParameter.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"

#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataTable.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"
#include "Offline/Mu2eUtilities/inc/ConversionSpectrum.hh"

namespace mu2e {

  //================================================================
  class LeadingLog : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),Comment("A SimParticleCollection with input stopped muons.")};
        fhicl::Atom<std::string> stoppingTargetMaterial{
        Name("stoppingTargetMaterial"),Comment("Material determines endpoint energy and muon life time.  Material must be known to the GlobalConstantsService."),"Al" };
        fhicl::DelegatedParameter captureProducts{Name("captureProducts"), Comment("A sequence of ParticleGenerator tools implementing capture products.")};
        fhicl::Atom<unsigned> verbosity{Name("verbosity"),0};
        fhicl::Atom<int> pdgId{Name("pdgId"),Comment("pdg id of daughter particle")};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit LeadingLog(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    void addParticles(StageParticleCollection* output, art::Ptr<SimParticle> mustop, double time, ParticleGeneratorTool* gen);
    //----------------------------------------------------------------
  private:
    const PDGCode::type electronId_ = PDGCode::e_minus; // for mass only
    double electronMass_;
    double endPointEnergy_;
    double endPointMomentum_;
    double muonLifeTime_;

    art::ProductToken<SimParticleCollection> const simsToken_;

    unsigned verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;
    RandomUnitSphere   randomUnitSphere_;
    ProcessCode process;
    int pdgId_;
    PDGCode::type pid;
    
    std::vector<std::unique_ptr<ParticleGeneratorTool>> muonCaptureGenerators_;
  };

  //================================================================
  LeadingLog::LeadingLog(const Parameters& conf)
    : EDProducer{conf}
    , electronMass_(GlobalConstantsHandle<ParticleDataTable>()->particle(electronId_).ref().mass().value())
    , endPointEnergy_()
    , endPointMomentum_ ()
    , muonLifeTime_{GlobalConstantsHandle<PhysicsParams>()->getDecayTime(conf().stoppingTargetMaterial())}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
    , randomUnitSphere_{eng_}
    , pdgId_(conf().pdgId())
  {
    produces<mu2e::StageParticleCollection>();
    pid = static_cast<PDGCode::type>(pdgId_);
    
    if (pid == PDGCode::e_minus) { 
      process = ProcessCode::mu2eCeMinusLeadingLog; 
      endPointEnergy_ = GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy(conf().stoppingTargetMaterial());
    } else if (pid == PDGCode::e_plus) { 
      process = ProcessCode::mu2eCeMinusLeadingLog; //TODO - new process code?
      endPointEnergy_ = GlobalConstantsHandle<PhysicsParams>()->getePlusEndpointEnergy(conf().stoppingTargetMaterial());
    }
    else {
      throw   cet::exception("BADINPUT")
        <<"LeadingLogGenerator::produce(): No process associated with chosen PDG id\n";
    }
    //endPointMomentum_ = endPointEnergy_*sqrt(1 - std::pow(electronMass_/endPointEnergy_,2));
    
    const auto capture_psets = conf().captureProducts.get<std::vector<fhicl::ParameterSet>>();
    for (const auto& pset : capture_psets) {
      muonCaptureGenerators_.push_back(art::make_tool<ParticleGeneratorTool>(pset));
      muonCaptureGenerators_.back()->finishInitialization(eng_, conf().stoppingTargetMaterial());
    }
  }

  //================================================================
  void LeadingLog::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};
    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto mus = stoppedMuMinusList(simh);
    

    for(const auto& mustop: mus) {

      const double time = mustop->endGlobalTime() + randExp_.fire(muonLifeTime_);
      //double rand = randFlat_.fire();

      
      for (const auto& gen : muonCaptureGenerators_) {
        addParticles(output.get(), mustop, time, gen.get());
      }


      if(mus.empty()) {
      throw   cet::exception("BADINPUT")
      <<"LeadingLog::produce(): no suitable stopped muon in the input SimParticleCollection\n";

      }

    
    }
  }
  
  //================================================================
  void LeadingLog::addParticles(StageParticleCollection* output,
                            art::Ptr<SimParticle> mustop,
                            double time,
                            ParticleGeneratorTool* gen)
  {
    auto daughters = gen->generate(); // this is were the energy should be used
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

DEFINE_ART_MODULE(mu2e::LeadingLog);
