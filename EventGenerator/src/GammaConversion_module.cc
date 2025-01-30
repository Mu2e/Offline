// Generates e+/e- from converted photons
// Uses the G4BetheHeitlerModel written into Mu2eUtilities::GammaConversionSpectrum
// This module throws an exception if no suitable photon is found.
//
// Expected workflow: Photons are generated and conversion events are filtered (e.g. with GammaConversionFilter + StoppedParticlesFinder), then resampled with this module
//
// M. MacKenzie (2024), based on GammaConvFlat_module.cc

#include <iostream>
#include <string>
#include <cmath>
#include <memory>

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandFlat.h"

#include "fhiclcpp/types/Atom.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/Mu2eUtilities/inc/GammaPairConversionSpectrum.hh"

// ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

namespace mu2e {

  //================================================================
  class GammaConversion : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),Comment("A SimParticleCollection with input converted photons.")};
      fhicl::Atom<int> material{Name("material"), Comment("Conversion material Z"), 13};
      fhicl::Atom<bool> correlateAnglesOverKE{Name("correlateAnglesOverKE"), Comment("Option to correlate angles over KE"), true};
      fhicl::Atom<int> verbosity{Name("verbosity"), Comment("Verbosity level"), 0};
      fhicl::Atom<bool> makeHistograms{Name("makeHistograms"), Comment("Make histograms of the conversion kinematics"), false};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit GammaConversion(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:

    art::ProductToken<SimParticleCollection> const simsToken_;
    const int material_;
    const int verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randFlat_;
    GammaPairConversionSpectrum spectrum_;
    GammaPairConversionSpectrum::elementData element_;
    const ProcessCode process_;
    const bool makeHistograms_;

    // histograms of the conversion kinematics
    TH1* hmomentum_;
    TH1* hCosz_;
    TH1* hEnergyElectron_;
    TH1* hEnergyPositron_;
    TH1* hMass_;
    TH2* hMassVsE_;
    TH1* hMassOverE_;
    TH1* hy_; // splitting function
  };

  //================================================================
  GammaConversion::GammaConversion(const Parameters& conf)
    : EDProducer{conf}
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , material_{conf().material()}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randFlat_{eng_}
    , spectrum_(&randFlat_, conf().correlateAnglesOverKE())
    , element_(spectrum_.GetElementMap()[material_])
    , process_(ProcessCode::mu2eGammaConversion)
    , makeHistograms_(conf().makeHistograms())
  {
    produces<mu2e::StageParticleCollection>();

    // make histograms of the conversion kinematics if requested
    if(makeHistograms_) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir("GammaConversion");
      hmomentum_       = tfdir.make<TH1F>("hmomentum", "Input photon momentum", 60,  0.,  120.  );
      hCosz_           = tfdir.make<TH1F>("hCosz", "Input photon cos(#theta_{z})", 200,  -1.,  1.  );
      hEnergyElectron_ = tfdir.make<TH1F>("hEnergyElectron", "Produced conversion electron energy", 60,  0.,  120.  );
      hEnergyPositron_ = tfdir.make<TH1F>("hEnergyPositron", "Produced conversion positron energy", 60,  0.,  120.  );
      hMass_           = tfdir.make<TH1F>("hMass"          , "M(e+e-) "           , 120,0.,120.);
      hMassVsE_        = tfdir.make<TH2F>("hMassVsE"       , "M(e+e-) vs. E;E (MeV); M(e+e-) (MeV/c^{2})"  , 120,0.,120.,120,0,120);
      hMassOverE_      = tfdir.make<TH1F>("hMassOverE"     , "M(e+e-)/E"          , 100, 0.,1);
      hy_              = tfdir.make<TH1F>("hy"             , "y = (ee-ep)/|pe+pp|", 100,-1.,1.);
    }
  }

  //================================================================
  void GammaConversion::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    // Sim particle collection from a generated photon event with a photon conversion
    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto gammas = simParticleList(simh,PDGCode::gamma,ProcessCode::conv);

    // If no suitable photons, print the particle collection and throw an error
    if(gammas.empty() || verbosity_ > 2) {
      std::cout << "GammaConversion_module::" << __func__ << ": Printing the particle collection:\n";
      const auto& sims = *simh;
      for(auto i = sims.begin(); i != sims.end(); ++i) {
        const auto& inpart = i->second;
        std::cout << "SimParticle: id = " << inpart.id() << " pdg = " << inpart.pdgId() << " createCode = " << inpart.creationCode()
                  << " stopCode = " << inpart.stoppingCode() << " E = " << inpart.endMomentum().e() << std::endl;
      }
    }
    if(gammas.empty()) {
      throw cet::exception("BADINPUT")
        << "GammaConversion::produce(): no suitable converted gamma in the input SimParticleCollection\n";
    }

    // Retrieve the first relevant photon from the collection, ignoring those from Bremsstrahlung or annihilations
    unsigned gamma_index = gammas.size();
    for(unsigned index = 0; index < gammas.size(); ++index) {
      const auto gamma = gammas.at(index);
      if(gamma->creationCode() == ProcessCode::eBrem || gamma->creationCode() == ProcessCode::annihil) continue;
      if(verbosity_ > 2) printf("--> Selecting photon at index %u\n", index);
      gamma_index = index;
      break;
    }

    // Check if a good photon is found
    if(gamma_index >= gammas.size()) {
      if(verbosity_ > 0) printf("GammaConversion::%s: No suitable photon found to convert, passing an empty event\n", __func__);
      event.put(std::move(output));
    }

    // Process the photon conversion
    const auto gamma = gammas.at(gamma_index);
    const auto p_g   = gamma->endMomentum();
    const auto t     = gamma->endGlobalTime();

    // Get a random conversion event given the input photon
    CLHEP::HepLorentzVector p_em, p_ep;
    spectrum_.fire(p_g, element_, p_em, p_ep);

    // Add the particles to the output
    output->emplace_back(gamma,
                         process_,
                         PDGCode::e_minus,
                         gamma->endPosition(),
                         p_em,
                         t
                         );

    output->emplace_back(gamma,
                         process_,
                         PDGCode::e_plus,
                         gamma->endPosition(),
                         p_ep,
                         t
                         );

    // print the conversion event if requested
    if (verbosity_ > 1){
      std::cout << "p(gamma) = (" << p_g .vect().x() << "," << p_g .vect().y() << "," << p_g .vect().z() << "," << p_g .e() << ")" << std::endl;
      std::cout << "p(e-)    = (" << p_em.vect().x() << "," << p_em.vect().y() << "," << p_em.vect().z() << "," << p_em.e() << ")" << std::endl;
      std::cout << "p(e+)    = (" << p_ep.vect().x() << "," << p_ep.vect().y() << "," << p_ep.vect().z() << "," << p_ep.e() << ")" << std::endl;
    }

    // fill the conversion kinematics histograms if requested
    if(makeHistograms_) {
      const double energy = p_g.e();
      const double mass = (p_em+p_ep).m();
      hmomentum_->Fill(energy);
      hCosz_->Fill((p_em+p_ep).cosTheta());
      hEnergyElectron_->Fill(p_em.e());
      hEnergyPositron_->Fill(p_ep.e());

      hMass_->Fill(mass);
      hMassVsE_->Fill(energy,mass);
      hMassOverE_->Fill(mass/energy);

      const CLHEP::Hep3Vector p = p_em.vect()+p_ep.vect();
      const double y = (p_em.e()-p_ep.e())/p.mag();
      hy_->Fill(y);
    }
    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::GammaConversion)
