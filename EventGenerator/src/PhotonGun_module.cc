// Generates photons
// Original author: Claudia Alvarez-Garcia
// Adapted by: Pawel Plesniak

// stdlib includes
#include <cmath>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

namespace mu2e {
  class PhotonGun : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<double> x{ Name("x"), Comment("x position of generated photon")};
      fhicl::Atom<double> y{ Name("y"), Comment("y position of generated photon")};
      fhicl::Atom<double> z{ Name("z"), Comment("z position of generated photon")};
      fhicl::OptionalAtom<double> px{ Name("px"), Comment("x momentum of generated photon")};
      fhicl::OptionalAtom<double> py{ Name("py"), Comment("y momentum of generated photon")};
      fhicl::OptionalAtom<double> deltax{ Name("deltax"), Comment("Difference in x position of generated photon, use to calculate targeted position")};
      fhicl::OptionalAtom<double> deltay{ Name("deltay"), Comment("Difference in y position of generated photon, use to calculate targeted position")};
      fhicl::OptionalAtom<double> deltaz{ Name("deltaz"), Comment("Difference in z position of generated photon, use to calculate targeted position")};
      fhicl::Atom<double> E{ Name("E"), Comment("Energy of generated photon")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit PhotonGun(const Parameters& conf);
    virtual void produce(art::Event& event);
  private:
    double x = 0.0, y = 0.0, z = 0.0;
    double px = 0.0, py = 0.0, pz = 0.0;
    double deltax = 0.0, deltay = 0.0, deltaz = 0.0;
    double E = 0.0;
  };

  PhotonGun::PhotonGun(const Parameters& conf):
    art::EDProducer(conf),
    x(conf().x()),
    y(conf().y()),
    z(conf().z()),
    E(conf().E()) {
      produces<GenParticleCollection>();
      if (E < std::numeric_limits<double>::epsilon())
        throw cet::exception("RANGE") << "Energy must be greater than zero, exiting.";
      px = conf().px() ? *conf().px() : 0;
      py = conf().py() ? *conf().py() : 0;
      if ((px*px + py*py) > (E*E))
        throw cet::exception("RANGE") << "magnitude of px and py is greater than E, exiting.";
      pz = std::sqrt(E*E - px*px - py*py);
      deltax = conf().deltax() ? *conf().deltax() : 0;
      deltay = conf().deltay() ? *conf().deltay() : 0;
      deltaz = conf().deltaz() ? *conf().deltaz() : 0;
      if ((px > std::numeric_limits<double>::epsilon() || py > std::numeric_limits<double>::epsilon()) && (deltax > std::numeric_limits<double>::epsilon() || deltay > std::numeric_limits<double>::epsilon() || deltaz > std::numeric_limits<double>::epsilon()))
        throw cet::exception("RANGE") << "Cannot specify both momentum and delta position, exiting.";
    };

  void PhotonGun::produce(art::Event& event) {
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    const CLHEP::Hep3Vector pos(x, y, z);
    CLHEP::Hep3Vector p(px, py, pz);
    if (std::abs(deltax) > std::numeric_limits<double>::epsilon() ||
        std::abs(deltay) > std::numeric_limits<double>::epsilon() ||
        std::abs(deltaz) > std::numeric_limits<double>::epsilon()) {
      std::cout << "PhotonGun: Using deltax, deltay, deltaz to calculate momentum." << std::endl;
      const CLHEP::Hep3Vector delta(deltax, deltay, deltaz);
      p = delta.unit() * E;
    }
    CLHEP::HepLorentzVector mom(p, E);
    output->push_back(GenParticle(PDGCode::gamma, GenId::particleGun, pos, mom, 0.));
    event.put(std::move(output));
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::PhotonGun)
