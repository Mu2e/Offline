// Generates outgoing photon or electron/positron pair  for either external or internal RPC
// Pions are generated with infinite lifetime, their real lifetime (proper time) is stored as a weight
// Andrei Gaponenko, 2013
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
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "Offline/SeedService/inc/SeedService.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/Mu2eUtilities/inc/PionCaptureSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/SimParticleGetTau.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"

// ROOT includes
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"

namespace mu2e {

  //================================================================
  class RPCGun : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
        fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),Comment("A SimParticleCollection with input stopped pions")};
        fhicl::Atom<unsigned> verbosity{Name("verbosity"),Comment("Add additional printout"), 0};
        fhicl::Atom<std::string> RPCType{Name("RPCType"),Comment("a process code, should be either RPCInternal or RPCExternal") };
        fhicl::DelegatedParameter spectrum{Name("spectrum"), Comment("Parameters for BinnedSpectrum")};
        fhicl::Atom<bool> pionDecayOff{Name("pionDecayOff"),Comment("Assume pion decay was turned off, produce event weights"), true};
        fhicl::Atom<bool> doHistograms{Name("doHistograms"),Comment("Produce debug histograms"), false};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit RPCGun(const Parameters& conf);

    virtual void produce(art::Event& event) override;
    void addParticles(StageParticleCollection* output,art::Ptr<SimParticle> pistop);
    double MakeEventWeight(art::Ptr<SimParticle> p);
    //----------------------------------------------------------------
  private:
    art::ProductToken<SimParticleCollection> const simsToken_;
    unsigned verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandExponential randExp_;
    CLHEP::RandFlat     randomFlat_;
    std::string RPCType_;
    BinnedSpectrum    spectrum_;
    bool pionDecayOff_;
    bool doHistograms_;
    RandomUnitSphere   randomUnitSphere_;
    CLHEP::RandGeneral randSpectrum_;

    ProcessCode process_;
    PionCaptureSpectrum pionCaptureSpectrum_;

    TH1F* _hnstops {nullptr};
    TH1F* _hmomentum {nullptr};
    TH1F* _hElecMom {nullptr};
    TH1F* _hElecPx {nullptr};
    TH1F* _hElecPy {nullptr};
    TH1F* _hElecPz {nullptr};
    TH1F* _hPosiMom {nullptr};
    TH1F* _hPosiPx {nullptr};
    TH1F* _hPosiPy {nullptr};
    TH1F* _hPosiPz {nullptr};
    TH1F* _hTime {nullptr};
    TH1F* _hTimeWt {nullptr};
    TH2F* _hStopXY {nullptr};
    TH1F* _hStopZ {nullptr};
    TH1F* _hWeight {nullptr};
    TH1F* _hMee;
    TH2F* _hMeeVsE;
    TH1F* _hMeeOverE;                   // M(ee)/E(gamma)
    TH1F* _hy;                          // splitting function

  };

  //================================================================
  RPCGun::RPCGun(const Parameters& conf)
    : EDProducer{conf}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randExp_{eng_}
    , randomFlat_{eng_}
    , RPCType_{conf().RPCType()}
    , spectrum_{BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>())}
    , pionDecayOff_{conf().pionDecayOff()}
    , doHistograms_{conf().doHistograms()}
    , randomUnitSphere_ {eng_}
    , randSpectrum_       {eng_, spectrum_.getPDF(),static_cast<int>(spectrum_.getNbins())}
    , pionCaptureSpectrum_{&randomFlat_,&randomUnitSphere_}
  {
    produces<mu2e::StageParticleCollection>();
    produces<mu2e::EventWeight>();
    if(RPCType_ == "mu2eExternalRPC") process_ = ProcessCode::mu2eExternalRPC;
    else if(RPCType_ == "mu2eInternalRPC") process_ = ProcessCode::mu2eInternalRPC;
    else {
      throw   cet::exception("BADINPUT")
        <<"RPCGun::produce(): no such process, must be mu2eInternalRPC or mu2eExternalRPC";
    }
    if ( doHistograms_ ) {
        art::ServiceHandle<art::TFileService> tfs;
        art::TFileDirectory tfdir = tfs->mkdir( "RPCGun" );

        _hnstops       = tfdir.make<TH1F>("hNStops", "N(pion stops) / event", 10,  -0.5,  9.5);
        _hmomentum     = tfdir.make<TH1F>("hmomentum", "Produced photon momentum", 100,  40.,  140.  );
        _hTime         = tfdir.make<TH1F>("hTime" , "pion time",100,0,1700);
        _hTimeWt       = tfdir.make<TH1F>("hTimeWt", "pion time, weighted",100,0,1700);
        _hStopXY       = tfdir.make<TH2F>("hStopXY", "Pion stop position;x (mm); y(mm)",100, -200, 200, 100, -200, 200);
        _hStopZ        = tfdir.make<TH1F>("hStopZ" , "Pion stop position;z (mm);", 200, 5000, 7000);
        _hWeight       = tfdir.make<TH1F>("hWeight" , "log10(weight)", 200, -10, 1.);
        if(RPCType_ == "mu2eInternalRPC"){
          _hElecMom  = tfdir.make<TH1F>("hElecMom" , "Produced electron momentum", 140,  0. , 140.);
          _hElecPx   = tfdir.make<TH1F>("hElecPx" , "Produced electron momentum Px", 140,  -140. , 140.);
          _hElecPy   = tfdir.make<TH1F>("hElecPy" , "Produced electron momentum Py", 140,  -140. , 140.);
          _hElecPz   = tfdir.make<TH1F>("hElecPz" , "Produced electron momentum Py", 140,  -140. , 140.);
          _hPosiMom  = tfdir.make<TH1F>("hPosiMom" , "Produced positron momentum", 140,  0. , 140.);
          _hPosiPx   = tfdir.make<TH1F>("hPosiPx" , "Produced positron momentum Px", 140,  -140. , 140.);
          _hPosiPy   = tfdir.make<TH1F>("hPosiPy" , "Produced positron momentum Py", 140,  -140. , 140.);
          _hPosiPz   = tfdir.make<TH1F>("hPosiPz" , "Produced positron momentum Pz", 140,  -140. , 140.);
          _hMee      = tfdir.make<TH1F>("hMee"     , "M(e+e-) "           , 200,0.,200.);
          _hMeeVsE   = tfdir.make<TH2F>("hMeeVsE"  , "M(e+e-) vs E"       , 200,0.,200.,200,0,200);
          _hMeeOverE = tfdir.make<TH1F>("hMeeOverE", "M(e+e-)/E "         , 200, 0.,1);
          _hy        = tfdir.make<TH1F>("hy"       , "y = (ee-ep)/|pe+pp|", 200,-1.,1.);
        }
      }
  }

  double RPCGun::MakeEventWeight(art::Ptr<SimParticle> part){
    if(verbosity_ > 1) printf("[RPCGun::%s] Evaluating the event weight from the particle proper time\n", __func__);
    const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>();
    const std::vector<int> decayOffCodes = {PDGCode::pi_plus, PDGCode::pi_minus};
    const double tau = SimParticleGetTau::calculate(part, decayOffCodes, gc);
    const double weight = exp(-tau);
    if(verbosity_ > 1) printf(" Tau = %.3f, Weight = %.3g, t(end) = %.3g\n", tau, weight, part->endGlobalTime());
    return weight;
  }

  //================================================================
  void RPCGun::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};
    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto pis = stoppedPiMinusList(simh);
    if(pis.empty()) {
      if(verbosity_ > 0) {
        printf("!!! RPCGun::%s: No pion stops found! Printing the sim collection:\n", __func__);
        const auto& sims = *simh;
        for(auto sim = sims.begin(); sim != sims.end(); ++sim) {
          std::cout << " " << sim->first << ":"
                    << " pdg = " << sim->second.pdgId()
                    << " code = " << sim->second.stoppingCode()
                    << " t = " << sim->second.endGlobalTime() << std::endl;
        }
      }
      throw   cet::exception("BADINPUT")
        <<"RPCGun::produce(): no suitable stopped pion in the input SimParticleCollection\n";
    }
    unsigned int randIn = randomFlat_.fireInt(pis.size());
    double time_weight = 1.0;
    if(pionDecayOff_) time_weight = MakeEventWeight(pis[randIn]);
    std::unique_ptr<EventWeight> pw(new EventWeight(time_weight));
    event.put(std::move(pw));
    addParticles(output.get(), pis[randIn]);
    event.put(std::move(output));
    if(doHistograms_) {
      _hnstops->Fill(pis.size());
    }
  }

  void RPCGun::addParticles(StageParticleCollection* output,
                            art::Ptr<SimParticle> pistop)
  {
    //Photon energy and four mom:
    double energy = spectrum_.sample(randSpectrum_.fire());
    const CLHEP::Hep3Vector p3 = randomUnitSphere_.fire(energy);
    const CLHEP::HepLorentzVector fourmom(p3, energy);
    if(verbosity_ > 0) printf("[RPCGun::%s] p_gamma = (%7.2f, %7.2f, %7.2f), t_gamma = %6.1f, pos = (%7.1f, %7.1f, %8.1f)\n", __func__,
                              p3.x(), p3.y(), p3.z(), pistop->endGlobalTime(),
                              pistop->endPosition().x(), pistop->endPosition().y(), pistop->endPosition().z());
    if(process_ == ProcessCode::mu2eExternalRPC){
     output->emplace_back(pistop,
                         process_,
                         PDGCode::gamma,
                         pistop->endPosition(),
                         fourmom,
                         pistop->endGlobalTime()
                         );

    } else if(process_ == ProcessCode::mu2eInternalRPC) {
      //Need to compute e-e+ pair momentum spectrum from the photon (use Kroll-Wada)
      CLHEP::HepLorentzVector mome, momp;
      pionCaptureSpectrum_.getElecPosiVectors(energy,mome,momp);
      output->emplace_back(pistop,
                           process_,
                           PDGCode::e_minus,
                           pistop->endPosition(),
                           mome,
                           pistop->endGlobalTime()
                           );

       output->emplace_back(pistop,
                           process_,
                           PDGCode::e_plus,
                           pistop->endPosition(),
                           momp,
                           pistop->endGlobalTime()
                           );

        if(doHistograms_){
          _hElecMom ->Fill(mome.vect().mag());
          _hElecPx ->Fill(mome.vect().x());
          _hElecPy ->Fill(mome.vect().y());
          _hElecPz ->Fill(mome.vect().z());
          _hPosiMom ->Fill(momp.vect().mag());
          _hPosiPx ->Fill(momp.vect().x());
          _hPosiPy ->Fill(momp.vect().y());
          _hPosiPz ->Fill(momp.vect().z());
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
          <<"RPCGun::produce(): no suitable process id\n";
      }
      if(doHistograms_) {
        const float weight = (pionDecayOff_) ? MakeEventWeight(pistop) : 1.;
        _hmomentum->Fill(energy);
        _hTime->Fill(pistop->endGlobalTime());
        _hTimeWt->Fill(pistop->endGlobalTime(), weight);
        _hStopXY->Fill(pistop->endPosition().x()+3904., pistop->endPosition().y());
        _hStopZ->Fill(pistop->endPosition().z());
        _hWeight->Fill((weight > 0.) ? std::log10(weight) : -1.e3);
      }

   }


  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::RPCGun)
