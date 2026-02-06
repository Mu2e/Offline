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
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"

#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"
#include "Offline/Mu2eUtilities/inc/ConversionSpectrum.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"

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
        fhicl::Atom<unsigned> verbosity{Name("verbosity"),0};
        fhicl::Atom<int> pdgId{Name("pdgId"),Comment("pdg id of daughter particle")};
        fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("Parameters for BinnedSpectrum)")};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit LeadingLog(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    void addParticles(StageParticleCollection* output, art::Ptr<SimParticle> mustop, double time);
    //----------------------------------------------------------------
  private:
    const PDGCode::type electronId_ = PDGCode::e_minus; // for mass only

    double muonLifeTime_;
    art::ProductToken<SimParticleCollection> const simsToken_;
    unsigned verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;
    ProcessCode process_;
    int pdgId_;
    PDGCode::type pid_;

    double _mass;
    BinnedSpectrum    spectrum_;
    std::unique_ptr<RandomUnitSphere>   randomUnitSphere_;
    std::unique_ptr<CLHEP::RandGeneral> randSpectrum_;
    double endPointEnergy_;

  };

  //================================================================
  LeadingLog::LeadingLog(const Parameters& conf)
    : EDProducer{conf}
    , muonLifeTime_{GlobalConstantsHandle<PhysicsParams>()->getDecayTime(conf().stoppingTargetMaterial())}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
    , pdgId_(conf().pdgId())
    ,_mass(GlobalConstantsHandle<ParticleDataList>()->particle(electronId_).mass())
    , spectrum_(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>()))
  {
    produces<mu2e::StageParticleCollection>();
    pid_ = static_cast<PDGCode::type>(pdgId_);

    if (pid_ == PDGCode::e_minus) {
      process_ = ProcessCode::mu2eCeMinusLeadingLog;
    } else if (pid_ == PDGCode::e_plus) {
      process_ = ProcessCode::mu2eCePlusLeadingLog;
    }
    else {
      throw   cet::exception("BADINPUT")
        <<"LeadingLogGenerator::produce(): No process associated with chosen PDG id : "<<pid_<<std::endl;
    }

    randomUnitSphere_ = std::make_unique<RandomUnitSphere>(eng_);
    randSpectrum_ = std::make_unique<CLHEP::RandGeneral>(eng_, spectrum_.getPDF(), spectrum_.getNbins());
  }

  //================================================================
  void LeadingLog::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};
    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto mus = stoppedMuMinusList(simh);

    if(mus.empty()) {
        throw   cet::exception("BADINPUT")
        <<"LeadingLog::produce(): no suitable stopped muon in the input SimParticleCollection\n";
    }

    for(const auto& mustop: mus) {
      const double time = mustop->endGlobalTime() + randExp_.fire(muonLifeTime_);
      addParticles(output.get(), mustop, time);
    }
    event.put(std::move(output));
  }

  //================================================================
  void LeadingLog::addParticles(StageParticleCollection* output,
                            art::Ptr<SimParticle> mustop,
                            double time)
  {

    double energy = spectrum_.sample(randSpectrum_->fire());

    const double p = sqrt((energy + _mass) * (energy - _mass));
    CLHEP::Hep3Vector p3 = randomUnitSphere_->fire(p);
    CLHEP::HepLorentzVector fourmom(p3, energy);

    output->emplace_back(mustop,
                       process_,
                       pid_,
                       mustop->endPosition(),
                       fourmom,
                       time
                       );

  }


  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::LeadingLog)
