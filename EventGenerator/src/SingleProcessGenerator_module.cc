// This module will be used to make electrons or positrons from any process originating from mu- stopped in any target material
//
// Sophie Middleton, 2021


#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandExponential.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
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
#include "Offline/GeometryService/inc/DetectorSystem.hh"

#include "KinKal/Trajectory/LoopHelix.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "TTree.h"
namespace mu2e {
  //================================================================
  class SingleProcessGenerator : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),
          Comment("A SimParticleCollection with muons stopping in the target.")};
      fhicl::Atom<std::string> stoppingTargetMaterial{
        Name("stoppingTargetMaterial"),
          Comment("Material determines muon life time, capture fraction, and particle spectra.\n"
                  "Only aluminum (Al) is supported, emisson spectra for other materials are not implemented.\n"),"Al" };
      fhicl::Atom<int> pdgId{Name("pdgId"),Comment("pdg id of daughter particle"),PDGCode::e_minus};
      fhicl::DelegatedParameter decayProducts{Name("decayProducts"), Comment("spectrum (and variables) to be generated")};
      fhicl::Atom<unsigned> verbosity{Name("verbosity"),0};
      fhicl::Atom<double> time_min{Name("time_min"),350};
      fhicl::Atom<double> time_max{Name("time_max"),1700};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit SingleProcessGenerator(const Parameters& conf);
    virtual void beginJob();
    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:
    double muonLifeTime_;
    art::ProductToken<SimParticleCollection> const simsToken_;
    unsigned verbosity_;
    double time_min_;
    double time_max_;
    int pdgId_;


    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;

    std::unique_ptr<ParticleGeneratorTool> Generator_;

    void addParticles(StageParticleCollection* output, art::Ptr<SimParticle> mustop, double time, ParticleGeneratorTool* gen);
    TTree* genTree;
    Float_t _maxr;
    Float_t _momT;
    Float_t _posT;
    Float_t _cosTheta;
    Float_t _time;
  };

  void SingleProcessGenerator::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    genTree  = tfs->make<TTree>("GenAna", "GenAna");
    genTree->Branch("maxr", &_maxr, "maxr/F");   
    genTree->Branch("momT", &_momT, "momT/F"); 
    genTree->Branch("posT", &_posT, "posT/F");
    genTree->Branch("cosTheta", &_cosTheta, "cosTheta/F");
    genTree->Branch("time", &_time, "time/F");
  }
  
  //================================================================
  SingleProcessGenerator::SingleProcessGenerator(const Parameters& conf)
    : EDProducer{conf}
    , muonLifeTime_{GlobalConstantsHandle<PhysicsParams>()->getDecayTime(conf().stoppingTargetMaterial())}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , time_min_{conf().time_min()}
    , time_max_{conf().time_max()}
    , pdgId_{conf().pdgId()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
  {
    produces<mu2e::StageParticleCollection>();
    if(verbosity_ > 0) {
      mf::LogInfo log("SingleProcessGenerator");
      log<<"stoppingTargetMaterial = "<<conf().stoppingTargetMaterial()
         <<", muon lifetime = "<<muonLifeTime_
         <<std::endl;
    }

    const auto pset = conf().decayProducts.get<fhicl::ParameterSet>();
    Generator_ = art::make_tool<ParticleGeneratorTool>(pset);
    Generator_->finishInitialization(eng_, conf().stoppingTargetMaterial());

    if(pdgId_==PDGCode::e_plus) {
      muonLifeTime_=0; //decay time already included for stopped muon(+) FIXME!!!
    }
    
  }

  //================================================================
  void SingleProcessGenerator::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto mus=(pdgId_==PDGCode::e_minus) ? stoppedMuMinusList(simh) : stoppedMuPlusList(simh);

     for(const auto& mustop: mus) {
      const double time = mustop->endGlobalTime() + randExp_.fire(muonLifeTime_);
      addParticles(output.get(), mustop, time, Generator_.get());
      }

    event.put(std::move(output));
  }

  //================================================================
  void SingleProcessGenerator::addParticles(StageParticleCollection* output,
                            art::Ptr<SimParticle> mustop,
                            double time,
                            ParticleGeneratorTool* gen)
  {
    auto daughters = gen->generate();
    for(const auto& d: daughters) {
     
      if(time > time_min_ and time < time_max_){

        //  make momentum and position vectors
        GeomHandle<DetectorSystem> det;
        ROOT::Math::XYZVectorF pos = XYZVectorF(det->toDetector(mustop->endPosition()));
        ROOT::Math::XYZTVector pos0(pos.x(), pos.y(), pos.z(), time); 
        ROOT::Math::PxPyPzMVector mom0(d.fourmom.x(), d.fourmom.y(), d.fourmom.z(), d.fourmom.t());

        // extract charge
        int charge = -1;
        if(pdgId_==PDGCode::e_plus) charge = 1;
        
        // extact field
        GeomHandle<BFieldManager> bfmgr;
        XYZVectorF pos3Vec = XYZVectorF(mustop->endPosition().x(),mustop->endPosition().y(),mustop->endPosition().z());
        ROOT::Math::XYZVector bnom(bfmgr->getBField(pos3Vec).x(),bfmgr->getBField(pos3Vec).y(),bfmgr->getBField(pos3Vec).z());

        // make the loophelix
        KinKal::LoopHelix lh(pos0, mom0, charge, bnom);
        // calculate rmax and add maxr to siminfo
        _maxr =sqrt(lh.cx()*lh.cx()+lh.cy()*lh.cy())+fabs(lh.rad());
        
        _momT =  sqrt(mom0.x()*mom0.x() + mom0.y()*mom0.y());
        _posT = sqrt(pos.x()*pos.x() + pos.y()*pos.y());
        _cosTheta = cos(atan2(_momT,mom0.z()));
        _time = time;
     
        output->emplace_back(mustop,
                             d.creationCode,
                             d.pdgId,
                             mustop->endPosition(),
                             d.fourmom,
                             time
                             );
      }
      genTree->Fill();
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SingleProcessGenerator)
