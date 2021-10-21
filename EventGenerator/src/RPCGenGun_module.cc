// Generates outgoing electron/positron pair from RPC
//
//  Sophie Middleton,2021

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

#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataTable.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "fhiclcpp/types/DelegatedParameter.h"

namespace mu2e {

  //================================================================
  class RPCGenGun : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),Comment("A SimParticleCollection with input stopped muons.")};
        fhicl::Atom<std::string> stoppingTargetMaterial{
        Name("stoppingTargetMaterial"),Comment("Material determines endpoint energy and pion life time.  Material pist be known to the GlobalConstantsService."),"Al" };
        fhicl::Atom<unsigned> verbosity{Name("verbosity"),0};
        fhicl::Atom<std::string> RPCType{
        Name("RPCType"),Comment("a process code, should be either RPCInternal or RPCExternal") };
        fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("Parameters for BinnedSpectrum)")};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit RPCGenGun(const Parameters& conf);

    virtual void produce(art::Event& event) override;
    void addParticles(StageParticleCollection* output,art::Ptr<SimParticle> pistop, FixedTimeMap time);
    //----------------------------------------------------------------
  private:
    const PDGCode::type electronId_ = PDGCode::e_minus; // for mass only
    double electronMass_;
    double pionLifeTime_;

    art::ProductToken<SimParticleCollection> const simsToken_;

    unsigned verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;
    RandomUnitSphere   randomUnitSphere_;
    ProcessCode process_;
    int pdgId_;
    PDGCode::type pid_;
    std::string RPCType_;
    BinnedSpectrum    spectrum_;
  };

  //================================================================
  RPCGenGun::RPCGenGun(const Parameters& conf)
    : EDProducer{conf}
    , electronMass_(GlobalConstantsHandle<ParticleDataTable>()->particle(electronId_).ref().mass().value())
    , pionLifeTime_{GlobalConstantsHandle<PhysicsParams>()->getDecayTime(conf().stoppingTargetMaterial())}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
    , randomUnitSphere_{eng_}
    , pdgId_(conf().pdgId())
    , RPCType_(conf().RPCType())
    , spectrum_(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>()))
  {
    produces<mu2e::StageParticleCollection>();
    produces<mu2e::EventWeight>();
    produces<mu2e::FixedTimeMap>();
    process_ = static_cast<PDGCode::type>(RPCType_);
  }

  //================================================================
  void RPCGenGun::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto pis = stoppedPiMinusList(simh);
    
    if(pis.empty()) {
      throw   cet::exception("BADINPUT")
        <<"RPCGenGun::produce(): no suitable stopped pion in the input SimParticleCollection\n";
    }
    
    
    std::unique_ptr<FixedTimeMap> timemap(new FixedTimeMap);

    ConditionsHandle<AcceleratorParams> accPar("ignored");
    double _mbtime = accPar->deBuncherPeriod;

    IO::StoppedParticleTauNormF stop;
    if (tmin_ > 0){
      while (true){
        const auto& tstop = stops_.fire();
        timemap->SetTime(protonPulse_->fire());
        if (tstop.t+timemap->time() < 0 || tstop.t+timemap->time() > tmin_){
          if (applySurvivalProbability_){
            double weight = exp(-tstop.tauNormalized)*survivalProbScaling_;
            if (weight > 1)
              std::cout << "WEIGHT TOO HIGH " << weight << " " << fmod(tstop.t,_mbtime) << std::endl;
            double rand = randomFlat_.fire();
            if (weight > rand){
              stop = tstop;
              break;
            }
          }else{
            stop = tstop;
            break;
          }
        }
      }
    }else{
      timemap->SetTime(protonPulse_->fire());
      if (applySurvivalProbability_){
        while (true){
          const auto& tstop = stops_.fire();
          double weight = exp(-tstop.tauNormalized)*survivalProbScaling_;
          double rand = randomFlat_.fire();
          if (weight > rand){
            stop = tstop;
            break;
          }
        }
      }else{
        const auto& tstop = stops_.fire();
        stop = tstop;
      }
    }
  //** USE the TOOL infrastructure **
    

    const auto pistop = pis.at(eng_.operator unsigned int() % pis.size());
    for(const auto& pistop: pis) {
      addParticles(output.get(), mustop, timemap);
    }
    event.put(std::move(output));
    event.put(std::move(timemap));
   
  }
  
  void RPCGenGun::addParticles(StageParticleCollection* output,
                            art::Ptr<SimParticle> pistop,
                            FixedTimeMap time)
  {
  
    if(process_ == ProcessCode::ExternalRPC){
     
     output->emplace_back(PDGCode::gamma,
                         process_, 
                         22,//FIXME
                         pistop->endPosition(),
                         CLHEP::HepLorentzVector{randopinitSphere_.fire(endPointMomentum_), endPointEnergy_},
                         pistop->endGlobalTime() + randExp_.fire(pionLifeTime_)
                         );

    } else if(process_ == ProcessCode::IntentalRPC) {
      output->emplace_back(PDGCode::e_minus,
                           , 
                           11,//FIXME
                           pistop->endPosition(),
                           CLHEP::HepLorentzVector{randopinitSphere_.fire(endPointMomentum_), endPointEnergy_},
                           pistop->endGlobalTime() + randExp_.fire(pionLifeTime_)
                           );
                           
       output->emplace_back(PDGCode::e_plus,
                           ProcessCode::InternalRPC, 
                           -11, //FIXME
                           pistop->endPosition(),
                           CLHEP::HepLorentzVector{randopinitSphere_.fire(endPointMomentum_), endPointEnergy_},
                           pistop->endGlobalTime() + randExp_.fire(pionLifeTime_)
                           );
                     event.put(std::move(output));
      
      } else {
        throw   cet::exception("BADINPUT")
        <<"RPCGenGun::produce(): no suitable process id\n";
      }
   }
  /*double RPCGun::generateEnergy() {
    return spectrum_.sample(randSpectrum_.fire());
  }*/

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::RPCGenGun);
