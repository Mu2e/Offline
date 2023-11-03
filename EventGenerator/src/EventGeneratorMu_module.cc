//-----------------------------------------------------------------------------
// keep is simple - generic generator module which uses plugins to
// generate various signal processes
//-----------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"

#include "canvas/Utilities/InputTag.h"

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
  class EventGeneratorMu : public art::EDProducer {
  public:
    struct Config {
      using Name   =fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> simpCollTag{Name("simpCollTag"  ),Comment("input SimParticleCollection")};
      fhicl::DelegatedParameter  generator  {Name("generator"    ),Comment("generator tool")};
      fhicl::Atom<int>           pdgCode    {Name("pdgCode"      ),Comment("PDG code"                   )};
      fhicl::Atom<double>        lifetime   {Name("lifetime"     ),Comment("lifetime"                   )};
      fhicl::Atom<unsigned>      verbosity  {Name("verbosity"    ),0};
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit EventGeneratorMu(const Parameters& conf);

    virtual void produce(art::Event& event) override;

//----------------------------------------------------------------
  private:
    double                                     lifetime_   ;
    int                                        pdgCode_    ;
    const art::InputTag                        simpCollTag_;
    unsigned                                   verbosity_  ;
    art::RandomNumberGenerator::base_engine_t& eng_        ;
    CLHEP::RandExponential                     randExp_    ;

    std::unique_ptr<ParticleGeneratorTool>     generator_  ;

    void addParticles(StageParticleCollection* output, art::Ptr<SimParticle> mustop, double time, ParticleGeneratorTool* gen);
  };

//================================================================
  EventGeneratorMu::EventGeneratorMu(const Parameters& conf) : EDProducer{conf}
    , lifetime_   {conf().lifetime   ()}
    , pdgCode_    {conf().pdgCode    ()}
    , simpCollTag_{conf().simpCollTag()}
    , verbosity_  {conf().verbosity  ()}
    , eng_        {createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_    {eng_}
  {
    consumes<SimParticleCollection>(simpCollTag_);
    produces<mu2e::StageParticleCollection>();

    const auto pset = conf().generator.get<fhicl::ParameterSet>();

    generator_      = art::make_tool<ParticleGeneratorTool>(pset);
    generator_->finishInitialization(eng_,"who was the genius to introduce this interface? - Let me know, thanks, Pasha");

    if(verbosity_ > 0) {
      mf::LogInfo log("EventGeneratorMu");
      log << ", lifetime = " << lifetime_ << std::endl;
    }
  }

  //================================================================
  void EventGeneratorMu::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simpch = event.getValidHandle<SimParticleCollection>(simpCollTag_);

    std::vector<art::Ptr<SimParticle>> list;

    if      (pdgCode_ ==   13) list = stoppedMuMinusList(simpch);
    else if (pdgCode_ ==  -13) list = stoppedMuPlusList (simpch);
    else if (pdgCode_ == -211) list = stoppedPiMinusList(simpch);

    for(const auto& simp: list) {
      double time = simp->endGlobalTime() + randExp_.fire(lifetime_);
      addParticles(output.get(), simp, time, generator_.get());
    }

    if(verbosity_ >= 9) {
      std::cout<<"EventGeneratorMu output: "<<*output<<std::endl;
    }

    event.put(std::move(output));
  }

  //================================================================
  void EventGeneratorMu::addParticles(StageParticleCollection* output,
                                      art::Ptr<SimParticle>    mustop,
                                      double                   time,
                                      ParticleGeneratorTool*    gen)
  {
    auto daughters = gen->generate();
    for(const auto& d: daughters) {
      output->emplace_back(mustop               ,
                           gen->processCode()   ,
                           d.pdgId              ,
                           mustop->endPosition(),
                           d.fourmom            ,
                           time                 );
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EventGeneratorMu)
