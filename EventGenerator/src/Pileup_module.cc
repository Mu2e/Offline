// This module implements Andy Edmonds' idea of generating the proper
// mix of "pileup" particles from muminus target stops in a single module.
// It is the MuStopProductsGun module re-cast to work with MC-history
// preserving StageParticles instead of GenParticles.
//
// Andrei Gaponenko, 2021

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

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataTable.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

#include <TH1F.h>
#include <TTree.h>
namespace mu2e {
  //================================================================
  class Pileup : public art::EDProducer {
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

      fhicl::DelegatedParameter captureProducts{Name("captureProducts"), Comment("A sequence of ParticleGenerator tools implementing capture products.")};
      fhicl::DelegatedParameter decayProducts{Name("decayProducts"), Comment("A sequence of ParticleGenerator tools implementing decay products.")};

      fhicl::Atom<unsigned> verbosity{Name("verbosity"),0};
      fhicl::Atom<bool> makeHistograms{Name("makeHistograms"),false};

    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit Pileup(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:
    Int_t nEv;
    Float_t pmag_gen;
    Float_t time_gen;
    Float_t x_gen;
    Float_t y_gen;
    Float_t z_gen;
    TTree*  _Ntup;
    
    double muonLifeTime_;
    double decayFraction_;

    art::ProductToken<SimParticleCollection> const simsToken_;

    unsigned verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randFlat_;
    CLHEP::RandExponential randExp_;
    bool makeHistograms_;
    std::vector<std::unique_ptr<ParticleGeneratorTool>> muonDecayGenerators_;
    std::vector<std::unique_ptr<ParticleGeneratorTool>> muonCaptureGenerators_;

    void addParticles(StageParticleCollection* output, art::Ptr<SimParticle> mustop, double time, ParticleGeneratorTool* gen);
  };

  //================================================================
  Pileup::Pileup(const Parameters& conf)
    : EDProducer{conf}
    , muonLifeTime_{GlobalConstantsHandle<PhysicsParams>()->getDecayTime(conf().stoppingTargetMaterial())}
    , decayFraction_{GlobalConstantsHandle<PhysicsParams>()->getDecayFraction(conf().stoppingTargetMaterial())}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randFlat_{eng_}
    , randExp_{eng_}
    , makeHistograms_{conf().makeHistograms()}
  {
    produces<mu2e::StageParticleCollection>();
    if(verbosity_ > 0) {
      mf::LogInfo log("Pileup");
      log<<"stoppingTargetMaterial = "<<conf().stoppingTargetMaterial()
         <<", muon lifetime = "<<muonLifeTime_
         <<", decay fraction = "<<decayFraction_
         <<std::endl;
    }
    if(makeHistograms_){
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "IPAGun" );

      _Ntup  = tfs->make<TTree>("GenTree", "GenTree");
	    _Ntup->Branch("nEv", &nEv , "nEv/I");	    
	    _Ntup->Branch("pmag_gen", &pmag_gen , "pmag_gen/F");
	    _Ntup->Branch("time_gen", &time_gen, "time_gen/F");
	    _Ntup->Branch("x_gen", &x_gen, "x_gen/F");
	    _Ntup->Branch("y_gen", &y_gen, "y_gen/F");
	    _Ntup->Branch("z_gen", &z_gen, "z_gen/F");
	  }
    if(conf().stoppingTargetMaterial() != "Al" and conf().stoppingTargetMaterial() != "IPA" ) {
      throw   cet::exception("NOT_IMPLEMENTED")
        <<"Pileup_module: emisson spectra for other than Al target are not impelmented\n";
    }

    const auto cap_psets = conf().captureProducts.get<std::vector<fhicl::ParameterSet>>();
    for (const auto& pset : cap_psets) {
      muonCaptureGenerators_.push_back(art::make_tool<ParticleGeneratorTool>(pset));
      muonCaptureGenerators_.back()->finishInitialization(eng_, conf().stoppingTargetMaterial());
    }

    const auto decay_psets = conf().decayProducts.get<std::vector<fhicl::ParameterSet>>();
    for (const auto& pset : decay_psets) {
      muonDecayGenerators_.push_back(art::make_tool<ParticleGeneratorTool>(pset));
      muonDecayGenerators_.back()->finishInitialization(eng_, conf().stoppingTargetMaterial());
    }
  }

  //================================================================
  void Pileup::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto mus = stoppedMuMinusList(simh);

    for(const auto& mustop: mus) {

      // decay or capture time for this muon, should
      // be the same for all its daughters
      const double time = mustop->endGlobalTime() + randExp_.fire(muonLifeTime_);
      std::cout<<mustop->endGlobalTime() <<" "<< randExp_.fire(muonLifeTime_)<<std::endl;
      double rand = randFlat_.fire();
      if (rand < decayFraction_) {
        for (const auto& gen : muonDecayGenerators_) {
          addParticles(output.get(), mustop, time, gen.get());
        }
      }
      else {
        for (const auto& gen : muonCaptureGenerators_) {
          addParticles(output.get(), mustop, time, gen.get());
        }
      }
      
    }

    
    if(verbosity_ >= 9) {
      std::cout<<"Pileup output: "<<*output<<std::endl;
      
    }

    event.put(std::move(output));
  }

  //================================================================
  void Pileup::addParticles(StageParticleCollection* output,
                            art::Ptr<SimParticle> mustop,
                            double time,
                            ParticleGeneratorTool* gen)
  {
    auto daughters = gen->generate();
    for(const auto& d: daughters) {
      //time = 1000;
      output->emplace_back(mustop,
                           d.creationCode,
                           d.pdgId,
                           mustop->endPosition(),
                           d.fourmom,
                           time
                           );
       if(makeHistograms_){
         double mag = sqrt(d.fourmom.px()*d.fourmom.px() + d.fourmom.py()*d.fourmom.py() + d.fourmom.pz()*d.fourmom.pz());
         pmag_gen = mag;
         time_gen = time;
         x_gen = mustop->endPosition().x();
         y_gen = mustop->endPosition().y();
         z_gen = mustop->endPosition().z();
         _Ntup->Fill();
      }
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::Pileup);
