// Generates e+/e- from converted photons
// Uses the G4BetheHeitlerModel written into Mu2eUtilities::GammaConversionSpectrum
// This module throws an exception if no suitable photon is found.
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

#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/Mu2eUtilities/inc/GammaPairConversionSpectrum.hh"

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
  {
    produces<mu2e::StageParticleCollection>();
  }

  //================================================================
  void GammaConversion::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto gammas = simParticleList(simh,PDGCode::gamma,ProcessCode::conv);

    // If no suitable photons, print the particle collection and throw an error
    if(gammas.empty()) {
      const auto& sims = *simh;
      std::cout << "GammaConversion_module::" << __func__ << ": Printing the particle collection:\n";
      for(auto i = sims.begin(); i != sims.end(); ++i) {
        const auto& inpart = i->second;
        std::cout << "SimParticle: pdg = " << inpart.pdgId() << " createCode = " << inpart.creationCode()
                  << " stopCode = " << inpart.stoppingCode() << std::endl;
      }
      throw cet::exception("BADINPUT")
        << "GammaConversion::produce(): no suitable converted gamma in the input SimParticleCollection\n";
    }

    // Retrieve a random photon from the collection
    const auto gamma = gammas.at(eng_.operator unsigned int() % gammas.size());
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

    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::GammaConversion)
