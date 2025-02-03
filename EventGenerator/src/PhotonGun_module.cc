// Generates photons
// Original author: Claudia Alvarez-Garcia
// Adapted by: Pawel Plesniak

// stdlib includes
#include <cmath>
#include <algorithm>
#include <iostream>

// exception handling
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/ParameterSet.h"

// Offline includes
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

using namespace std;
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
      fhicl::Atom<double> E{ Name("E"), Comment("Energy of generated photon")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit PhotonGun(const Parameters& conf);
    virtual void produce(art::Event& event);
  private:
    double x = 0.0, y = 0.0, z = 0.0;
    double px = 0.0, py = 0.0, pz = 0.0;
    double E = 0.0;
  };

  PhotonGun::PhotonGun(const Parameters& conf):
    art::EDProducer(conf),
    x(conf().x()),
    y(conf().y()),
    z(conf().z()),
    E(conf().E()) {
      produces<GenParticleCollection>();
      px = conf().px() ? *conf().px() : 0;
      px = conf().py() ? *conf().py() : 0;
      if ((px*px + py*py) > (E*E))
        throw cet::exception("RANGE") << "magnitude of px and py is greater than E, exiting.";
      pz = sqrt(E*E - px*px - py*py);
    };

  void PhotonGun::produce(art::Event& event) {
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    const CLHEP::Hep3Vector pos(x, y, z);
    const CLHEP::Hep3Vector p(px, py, pz);
    CLHEP::HepLorentzVector mom(p, E);
    output->push_back(GenParticle(PDGCode::gamma, GenId::particleGun, pos, mom, 0.));
    event.put(std::move(output));
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::PhotonGun)
