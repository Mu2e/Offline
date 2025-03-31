// Output stopped antiproton events from an input antiproton sim particle list
// Michael MacKenize, 2024


#include <iostream>
#include <string>
#include <cmath>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"

#include "TH1.h"

namespace mu2e {

  //================================================================
  class AntiprotonResampling : public art::EDProducer {
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"), Comment("Input sim particle collection")};
      fhicl::Atom<int> verbosity{Name("verbosity"), Comment("Verbosity level"), 0};
      fhicl::Atom<bool> makeHistograms{Name("makeHistograms"), Comment("Make histograms of the conversion kinematics"), false};
    };


    PDGCode::type       pdgId_;
    ProcessCode         processCode_;
    double              mass_;

    art::ProductToken<SimParticleCollection> const simsToken_;

    int                 verbosity_;

    //-----------------------------------------------------------------------------
    // histogramming
    //-----------------------------------------------------------------------------
    bool    makeHistograms_;
    TH1*   _hX;
    TH1*   _hY;
    TH1*   _hZ;
    TH1*   _hR;
    TH1*   _hTime;

  private:

  public:
    using Parameters= art::EDProducer::Table<Config>;
    explicit AntiprotonResampling(const Parameters& conf);

    virtual void produce(art::Event& event);
  };

  //================================================================
  AntiprotonResampling::AntiprotonResampling(const Parameters& conf)
    : EDProducer{conf}
    , pdgId_(PDGCodeDetail::anti_proton)
    , processCode_(ProcessCode::mu2eAntiproton)
    , mass_(GlobalConstantsHandle<ParticleDataList>()->particle(pdgId_).mass())
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_(conf().verbosity())
    , makeHistograms_(conf().makeHistograms())
  {
    produces<mu2e::StageParticleCollection>();

    if(verbosity_ > 0) {
      std::cout<<"AntiprotonResampling: using process code " << processCode_ << std::endl;
    }


    if(makeHistograms_) {
      art::ServiceHandle<art::TFileService> tfs;
      _hX      = tfs->make<TH1F>("hX"     , "X"           ,  100,     0.,  100.);
      _hY      = tfs->make<TH1F>("hY"     , "Y"           ,  100,     0.,  100.);
      _hZ      = tfs->make<TH1F>("hZ"     , "Z"           ,  500,  5400., 6400.);
      _hR      = tfs->make<TH1F>("hR"     , "R"           ,  100,     0.,  100.);
      _hTime   = tfs->make<TH1F>("hTime"  , "Time"        ,  400,     0., 2000.);
    }
  }

  //================================================================
  void AntiprotonResampling::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    // Input sim particle collection
    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    auto sims = simParticleList(simh, pdgId_);

    if(sims.empty()) {
      throw cet::exception("BADINPUT") << "AntiprotonResampling::" << __func__ << ": No input sim particles found\n";
      // if(verbosity_ > 0) printf("AntiprotonResampling::%s: No antiprotons found, inserting an empty particle list\n", __func__);
      // event.put(std::move(output));
    }

    // Only take the first sim particle (only one per event is expected)
    const auto sim = sims[0];

    // Annihilation position
    const CLHEP::Hep3Vector pos(sim->endPosition());

    // start at a stop
    const double energy = mass_;
    CLHEP::Hep3Vector p3(0., 0., 0.);
    CLHEP::HepLorentzVector fourmom(p3, energy);
    output->emplace_back(sim,
                         processCode_,
                         pdgId_,
                         pos,
                         fourmom,
                         sim->endGlobalTime());

    event.put(std::move(output));

    if(verbosity_ > 1) {
      printf("AntiprotonResampling::%s: Add antiproton with (x,y,z,t) = (%7.1f, %5.1f, %8.1f, %6.1f)\n",
             __func__, pos.x(), pos.y(), pos.z(), sim->endGlobalTime());
    }
    //-----------------------------------------------------------------------------
    // if requested, fill histograms
    //-----------------------------------------------------------------------------
    if(makeHistograms_) {
      _hX->Fill(pos.x() + 3904.);
      _hY->Fill(pos.y());
      _hZ->Fill(pos.z());
      _hTime->Fill(sim->endGlobalTime());
      const double r = std::sqrt(std::pow(pos.x() + 3904., 2) + std::pow(pos.y(), 2));
      _hR->Fill(r);
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::AntiprotonResampling)
