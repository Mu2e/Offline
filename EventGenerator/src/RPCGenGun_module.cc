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
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"
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
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/Mu2eUtilities/inc/PionCaptureSpectrum.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/Mu2eUtilities/inc/SimParticleGetTau.hh"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"

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
        fhicl::Atom<std::string> RPCType{Name("RPCType"),Comment("a process code, should be either RPCInternal or RPCExternal") };
        fhicl::Sequence<int> decayOffPDGCodes{Name("decayOffPDGCodes"),Comment("decayOffPDGCodes")};
        fhicl::Sequence<art::InputTag> hitCollections {Name("hitCollections"), Comment("A list of StepPointMCCollection-s")};
        fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("Parameters for BinnedSpectrum")};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit RPCGenGun(const Parameters& conf);

    virtual void produce(art::Event& event) override;
    void addParticles(StageParticleCollection* output,art::Ptr<SimParticle> pistop);
    void MakeEventWeight(art::Event& event);
    //----------------------------------------------------------------
  private:
    double pionLifeTime_;
    art::ProductToken<SimParticleCollection> const simsToken_;
    unsigned verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;
    CLHEP::RandFlat     randomFlat_;
    std::string RPCType_;
    std::vector<int> decayOffPDGCodes_;
    std::vector<art::InputTag> hitColls_;
    BinnedSpectrum    spectrum_;
    RandomUnitSphere   randomUnitSphere_;
    int nbins_;
    CLHEP::RandGeneral randSpectrum_;
    


    double            tmin_ = -1;
    double            czmin_ = -1;
    double            czmax_ = 1;
    double            phimin_ = 0;
    double            phimax_ = CLHEP::twopi;
    ProcessCode process_;
    PionCaptureSpectrum pionCaptureSpectrum_;
    
    TH1F* _hmomentum;
    TH1F* _hElecMom {nullptr};
    TH1F* _hPosiMom {nullptr};
    TH1F* _hMee;
    TH2F* _hMeeVsE;
    TH1F* _hMeeOverE;                   // M(ee)/E(gamma)
    TH1F* _hy;                          // splitting function
    
    bool doHistograms_;
  };

  //================================================================
  RPCGenGun::RPCGenGun(const Parameters& conf)
    : EDProducer{conf}
    , pionLifeTime_{GlobalConstantsHandle<PhysicsParams>()->getDecayTime(conf().stoppingTargetMaterial())}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
    , randomFlat_{eng_}
    , RPCType_{conf().RPCType()}
    , decayOffPDGCodes_{conf().decayOffPDGCodes()}
    , hitColls_{conf().hitCollections()}
    , spectrum_{BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>())} 
    , randomUnitSphere_ {eng_, czmin_,czmax_,phimin_,phimax_}
    , nbins_{static_cast<int>(spectrum_.getNbins())}
    , randSpectrum_       {eng_, spectrum_.getPDF(),nbins_}
    , pionCaptureSpectrum_{&randomFlat_,&randomUnitSphere_}
  {
    produces<mu2e::StageParticleCollection>();
    produces<mu2e::EventWeight>();
    if(RPCType_ == "ExternalRPC") process_ = ProcessCode::ExternalRPC;
    else if(RPCType_ == "InternalRPC") process_ = ProcessCode::InternalRPC;
    else {
      throw   cet::exception("BADINPUT")
        <<"RPCGenGun::produce(): no such process, must be InternalRPC or ExternalRPC";
    } //TODO - replace with static_cast<ProcessCode::type>(RPCType_);  - got errors  for some reason so simplified

    //randomUnitSphere_ = new RandomUnitSphere(eng_);//,czmin_,czmax_,phimin_,phimax_);
    //randSpectrum_ = new CLHEP::RandGeneral(eng_, spectrum_.getPDF(), spectrum_.getNbins());
    //pionCaptureSpectrum_ = new PionCaptureSpectrum(randomFlat_,randomUnitSphere_); //make pion capture spectrum
    if ( doHistograms_ ) {
        art::ServiceHandle<art::TFileService> tfs;
        art::TFileDirectory tfdir = tfs->mkdir( "RPCGun" );

        _hmomentum     = tfdir.make<TH1F>( "hmomentum", "Produced photon momentum", 100,  40.,  140.  );

        if(RPCType_ == "InternalRPC"){
          _hElecMom  = tfdir.make<TH1F>("hElecMom" , "Produced electron momentum", 140,  0. , 140.);
          _hPosiMom  = tfdir.make<TH1F>("hPosiMom" , "Produced positron momentum", 140,  0. , 140.);
          _hMee      = tfdir.make<TH1F>("hMee"     , "M(e+e-) "           , 200,0.,200.);
          _hMeeVsE   = tfdir.make<TH2F>("hMeeVsE"  , "M(e+e-) vs E"       , 200,0.,200.,200,0,200);
          _hMeeOverE = tfdir.make<TH1F>("hMeeOverE", "M(e+e-)/E "         , 200, 0.,1);
          _hy        = tfdir.make<TH1F>("hy"       , "y = (ee-ep)/|pe+pp|", 200,-1.,1.);
        }
      }
  }
  
  void RPCGenGun::MakeEventWeight(art::Event& event){ //TODO - make this a utility
    std::vector<StepPointMCCollection> spMCColls;
    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_); 
    for ( const auto& iColl : hitColls_ ){
      auto spColl = event.getValidHandle<StepPointMCCollection>(iColl); 
      spMCColls.push_back( *spColl );
    }
    double weight = 0.;

    for(const auto& p : *simh) {
      if(p.second.daughters().empty()) {
          art::Ptr<SimParticle> pp(simh, p.first.asUint());
          double sumTime = SimParticleGetTau::calculate(pp, spMCColls, decayOffPDGCodes_);
          weight += exp(-sumTime);
        }
    }
    std::unique_ptr<EventWeight> pw(new EventWeight(weight));
    event.put(std::move(pw)); 
  }

  //================================================================
  void RPCGenGun::produce(art::Event& event) {
    std::cout<<"[RPCGenGun:: produce] start "<<std::endl;
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_); 
    const auto pis = stoppedPiMinusList(simh);
    MakeEventWeight(event);
    if(pis.empty()) {
      throw   cet::exception("BADINPUT")
        <<"RPCGenGun::produce(): no suitable stopped pion in the input SimParticleCollection\n";
    }
    for(const auto& pistop: pis) {
      addParticles(output.get(), pistop);
    }
    event.put(std::move(output)); //TODO - delete once addParticle function works
    std::cout<<"[RPCGenGun:: produce] end "<<std::endl;
  }
  
  void RPCGenGun::addParticles(StageParticleCollection* output,
                            art::Ptr<SimParticle> pistop)
  {
    std::cout<<"[RPCGenGun:: add] start"<<std::endl;
    //Photon energy and four mom:
    double energy = spectrum_.sample(randSpectrum_.fire());
    CLHEP::Hep3Vector p3 = randomUnitSphere_.fire(energy);
    CLHEP::HepLorentzVector fourmom(p3, energy);
    if(process_ == ProcessCode::ExternalRPC){
     
     output->emplace_back(pistop,
                         process_, 
                         PDGCode::gamma,
                         pistop->endPosition(),
                         fourmom,
                         pistop->endGlobalTime() + randExp_.fire(pionLifeTime_) //TODO - need to add event weight
                         );

    } else if(process_ == ProcessCode::InternalRPC) {
      //Need to compute e-e+ pair momentum spectrum from the photon (use Kroll-Wada)
      CLHEP::HepLorentzVector mome, momp;
      pionCaptureSpectrum_.getElecPosiVectors(energy,mome,momp);
      output->emplace_back(pistop,
                           process_, 
                           PDGCode::e_minus,
                           pistop->endPosition(),
                           mome,
                           pistop->endGlobalTime() + randExp_.fire(pionLifeTime_)
                           );
                           
       output->emplace_back(pistop,
                           process_, 
                           PDGCode::e_plus, 
                           pistop->endPosition(),
                           momp,
                           pistop->endGlobalTime() + randExp_.fire(pionLifeTime_)
                           );
                  
        if(doHistograms_){
          _hElecMom ->Fill(mome.vect().mag());
          _hPosiMom ->Fill(momp.vect().mag());

          double mee = (mome+momp).m();
          _hMee->Fill(mee);
          _hMeeVsE->Fill(energy,mee);
          _hMeeOverE->Fill(mee/energy);

          CLHEP::Hep3Vector p = mome.vect()+momp.vect();
          double y = (mome.e()-momp.e())/p.mag();

          _hy->Fill(y);
        }
      } else {
        throw   cet::exception("BADINPUT")
        <<"RPCGenGun::produce(): no suitable process id\n";
      }
      _hmomentum->Fill(energy);
      std::cout<<"[RPCGenGun:: add] end"<<std::endl;
   }


  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::RPCGenGun);
