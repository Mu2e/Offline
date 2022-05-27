// This module will be used to make electrons from any process originating from mu- stopped in any target material
//
// Sophie Middleton, 2021

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "art_root_io/TFileService.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Vector/ThreeVector.h"
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
#include "Offline/DataProducts/inc/GenVector.hh"

namespace mu2e {
  //================================================================
  class MuplusProcessGenerator : public art::EDProducer {
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
      fhicl::DelegatedParameter decayProducts{Name("decayProducts"), Comment("spectrum (and variables) to be generated")};
      fhicl::Atom<unsigned> verbosity{Name("verbosity"),0};

    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit MuplusProcessGenerator(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:
    double muonLifeTime_;
    art::ProductToken<SimParticleCollection> const simsToken_;
    unsigned verbosity_;

    TTree* _Ntup;
    Int_t _genPdgId;
    Int_t _genCrCode;
    Float_t _genPz;
    Float_t _genPosZ;
    Float_t _genPosR;
    Float_t _genTime;
    Float_t _genP; 

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;
  
    std::unique_ptr<ParticleGeneratorTool> Generator_;

    void addParticles(StageParticleCollection* output, art::Ptr<SimParticle> mustop, double time, ParticleGeneratorTool* gen);
  };

  //================================================================
  MuplusProcessGenerator::MuplusProcessGenerator(const Parameters& conf)
    : EDProducer{conf}
  // , muonLifeTime_{GlobalConstantsHandle<PhysicsParams>()->getDecayTime(conf().stoppingTargetMaterial())}
    , muonLifeTime_{0} //for mu+ stops only
 
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
  {
    produces<mu2e::StageParticleCollection>();
    if(verbosity_ > 0) {
      mf::LogInfo log("MuplusProcessGenerator");
      log<<"stoppingTargetMaterial = "<<conf().stoppingTargetMaterial()
         <<", muon lifetime = "<<muonLifeTime_
         <<std::endl;
    }

    const auto pset = conf().decayProducts.get<fhicl::ParameterSet>();
    Generator_ = art::make_tool<ParticleGeneratorTool>(pset);
    Generator_->finishInitialization(eng_, conf().stoppingTargetMaterial());

   art::ServiceHandle<art::TFileService> tfs;
   _Ntup = tfs->make<TTree>("GenAna","GenAna");
   _Ntup -> Branch("genId",&_genPdgId,"genId/I");
  _Ntup -> Branch("genCrCode",&_genCrCode,"genCrCode/I");
 _Ntup -> Branch("genPz",&_genPz,"genPz/F");
 _Ntup -> Branch("genPosZ",&_genPosZ,"genPosZ/F");
 _Ntup -> Branch("genPosR",&_genPosR,"genPosR/F");
 _Ntup -> Branch("genTime",&_genTime,"genTime/F");
  _Ntup -> Branch("genP",&_genP,"genP/F");
    
  }

  //================================================================
  void MuplusProcessGenerator::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto mus = stoppedMuPlusList(simh);

    for(const auto& mustop: mus) {

      const double time = mustop->endGlobalTime() + randExp_.fire(muonLifeTime_);
      addParticles(output.get(), mustop, time, Generator_.get());

    }

    event.put(std::move(output));
  }

  //================================================================
  void MuplusProcessGenerator::addParticles(StageParticleCollection* output,
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

      _genPdgId=d.pdgId;
      _genCrCode=d.creationCode;
      _genPosZ= mustop->endPosition().z();
      _genPosR= sqrt((mustop->endPosition().x()+3904)*(mustop->endPosition().x()+3904)+(mustop->endPosition().y())*(mustop->endPosition().y()));
      _genTime=time;
      _genP=sqrt(d.fourmom.px()*d.fourmom.px()+d.fourmom.py()*d.fourmom.py()+d.fourmom.pz()*d.fourmom.pz());
      _genPz=d.fourmom.pz();
      //    printf("the z position is %f\n",_genPosZ);
      //   printf("the R position is %f\n",_genPosR);
      _Ntup->Fill();

    }

   

  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::MuplusProcessGenerator);
