// Generates e+/e- from converted photons.
// This FAST generator uses the following generation technique:
// 0) read converted photon momentum from SimParticle collection
// 1) e- momentum is generated with a flat spectrum from 0 to the photon momentum
// 2) e- direction with respect to the photon direction is taken fixed at (electron mass)/(gamma momentum)
// 3) e- momentum angle in the plane transerse to the photon momentum is uniform from 0 to 2pi
// 4) e- momentum wrt to photon momentum is then rotated to get it in the Mu2e frame
// 5) e+ momentum is obtained imposing momentum conservation
// 6) e- and e+ are then added from the input SimParticleCollection.
// See also doc-db 47028
// This module throws an exception if no suitable photon is found.
//
// S. Di Falco, 2023

#include <iostream>
#include <string>
#include <cmath>
#include <memory>

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "fhiclcpp/types/Atom.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"

namespace mu2e {

  //================================================================
  class GammaConvFlat : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> inputSimParticles{Name("inputSimParticles"),Comment("A SimParticleCollection with input converted photons.")};
      fhicl::Atom<unsigned> verbosity{Name("verbosity")};
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit GammaConvFlat(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:

    double particleMass_;

    art::ProductToken<SimParticleCollection> const simsToken_;

    unsigned verbosity_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randFlat_;
    ProcessCode process;
  };

  //================================================================
  GammaConvFlat::GammaConvFlat(const Parameters& conf)
    : EDProducer{conf}
    , particleMass_(GlobalConstantsHandle<ParticleDataList>()->particle(static_cast<PDGCode::type>(PDGCode::e_minus)).mass())
    , simsToken_{consumes<SimParticleCollection>(conf().inputSimParticles())}
    , verbosity_{conf().verbosity()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randFlat_{eng_}
  {
    produces<mu2e::StageParticleCollection>();
    process = ProcessCode::mu2eGammaConversion;
  }

  //================================================================
  void GammaConvFlat::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};

    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto gammas = simParticleList(simh,PDGCode::gamma,ProcessCode::conv);

    if(gammas.empty()) {
      throw   cet::exception("BADINPUT")
        <<"GammaConvFlat::produce(): no suitable converted gamma in the input SimParticleCollection\n";
    }

    // Take the converted gamma as SimParticle
    const auto gammapart = gammas.at(eng_.operator unsigned int() % gammas.size());
    const auto p4_gamma = gammapart->endMomentum();
    const auto p_gamma   = p4_gamma.vect().mag();
    const auto t_gammastop = gammapart->endGlobalTime();

    // Start generating the e-:
    // extract flat momentum between 0 and p_gamma
    double p_eminus = randFlat_.fire(0., p_gamma);
    double E_eminus = sqrt(particleMass_*particleMass_ + p_eminus*p_eminus);
    // extract flat between 0 and 2pi the phi angle wrt gamma
    double phi_eminus = randFlat_.fire(0., CLHEP::twopi);
    // Fix to Mass_electron/Egamma radians the theta angle wrt to gamma
    double theta_eminus = particleMass_/p_gamma;
    // rotate so that the photon direction is the new z axis
    const auto dir_gamma = p4_gamma.vect().unit();
    double px_eminus = p_eminus * sin(theta_eminus) * cos (phi_eminus);
    double py_eminus = p_eminus * sin(theta_eminus) * sin (phi_eminus);
    double pz_eminus = p_eminus * cos(theta_eminus);
    CLHEP::HepLorentzVector p4_eminus_temp(px_eminus, py_eminus, pz_eminus, E_eminus);
    CLHEP::HepLorentzVector p4_eminus= p4_eminus_temp.rotateUz(dir_gamma);

    output->emplace_back(gammapart,
                         process,
                         PDGCode::e_minus,
                         gammapart->endPosition(),
                         p4_eminus,
                         t_gammastop
                         );

    // Now impose momentum conservation to find positron momentum

    double px_eplus = -px_eminus;
    double py_eplus = -py_eminus;
    double pz_eplus = p_gamma-pz_eminus;
    double p_eplus  = sqrt( px_eplus*px_eplus +
                            py_eplus*py_eplus +
                            pz_eplus*pz_eplus);
    double E_eplus = sqrt(particleMass_*particleMass_ + p_eplus*p_eplus );
    CLHEP::HepLorentzVector p4_eplus_temp(px_eplus, py_eplus, pz_eplus, E_eplus);
    CLHEP::HepLorentzVector p4_eplus= p4_eplus_temp.rotateUz(dir_gamma);

    output->emplace_back(gammapart,
                         process,
                         PDGCode::e_plus,
                         gammapart->endPosition(),
                         p4_eplus,
                         t_gammastop
                         );

    if (verbosity_>1){
      std::cout << "p(gamma)=(" << p4_gamma.vect().x() << "," << p4_gamma.vect().y() << "," << p4_gamma.vect().z() << "," << p_gamma << ")" << std::endl;
      std::cout << "p(eminus)=(" << p4_eminus.vect().x() << "," << p4_eminus.vect().y() << "," << p4_eminus.vect().z() << "," << E_eminus << ")" << std::endl;
      std::cout << "p(eplus)=(" << p4_eplus.vect().x() << "," << p4_eplus.vect().y() << "," << p4_eplus.vect().z() << "," << E_eplus << ")" << std::endl;
    }

    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::GammaConvFlat)
