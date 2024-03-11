// This generator throws particles of the same type and the same
// kinematics in each event, i.e. a delta-function distribution in the
// phase space.  This is exactly what is needed for some simulations,
// in particular those using G4study.  The explicit module name will
// hopefully help to prevent future complexification of this code.
// This generator is so basic that it does not even need random numbers.
// It is fully configurable via fcl, without the need for an
// additional txt file.
//
// Andrei Gaponenko, 2024

#include <string>
#include <memory>
#include <vector>
#include <cmath>

#include "cetlib_except/exception.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/TupleAs.h"

#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Principal/Event.h"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"

namespace mu2e {

  class deltaFunctionGun : public art::SharedProducer {
    public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      using Hep3Vector_t = CLHEP::Hep3Vector(double, double, double);

      fhicl::Atom<std::string> particleName {
        Name("particleName"),
        Comment("Particle name as defined in the PDGCode enum-to-string")
      };

      fhicl::TupleAs<Hep3Vector_t> position {
        Name("position"),
        Comment("The point from which the particles are launched, in the GenParticle coordinate system.\n"
                "This is usually the same as the Mu2e coordinate system, but the latter is not defined\n"
                "for simplified geometries.  The units are mm."
                )
      };

      fhicl::TupleAs<Hep3Vector_t> direction {
        Name("direction"),
        Comment("A non-zero 3-vector ponting in the direction of the particle momentum, in the GenParticle\n"
                "coordinate system."
                )
      };

      fhicl::OptionalAtom<double> momentum {
        Name("momentum"),
        Comment("The magnitude of momentum of the generated particles, in MeV/c.\n"
                "This setting is mutially exclusive with kineticEnergy.")
      };

      fhicl::OptionalAtom<double> kineticEnergy {
        Name("kineticEnergy"),
        Comment("Kinetic energy of the generated particles, in MeV/c2.\n"
                "This setting is mutially exclusive with momentum.")
      };
    };

    using Parameters = art::SharedProducer::Table<Config>;
    explicit deltaFunctionGun(const Parameters& conf, const art::ProcessingFrame&);

    void produce(art::Event& event, const art::ProcessingFrame&) override;

  private:
    const GenParticle particle_;
    GenParticle define_particle(const Parameters& conf);
  };

  //================================================================
  deltaFunctionGun::deltaFunctionGun(const Parameters& conf, const art::ProcessingFrame&)
    : SharedProducer{conf}
    , particle_{define_particle(conf)}
    {
      async<art::InEvent>();
      produces<GenParticleCollection>();
    }

  //================================================================
  void deltaFunctionGun::produce(art::Event& event, const art::ProcessingFrame&) {
    using namespace std;
    auto output{make_unique<GenParticleCollection,initializer_list<GenParticle>>({particle_})};
    event.put(move(output));
  }

  //================================================================
  GenParticle deltaFunctionGun::define_particle(const Parameters& conf) {

    PDGCode pdg(conf().particleName());
    GlobalConstantsHandle<ParticleDataList> pdt;
    const double mass = pdt->particle(pdg.id()).mass();

    //----------------------------------------------------------------
    // figure out the scalar momentum and energy
    double pmag = 0., etot = 0.;
    if(conf().momentum()) {
      if(conf().kineticEnergy()) {
        throw cet::exception("BADCONFIG")
          << "deltaFunctionGun: ERROR: both momentum and kineticEnergy are specified"
          << "\n";
      }

      pmag = *conf().momentum();
      if(pmag < 0.) {
        throw cet::exception("BADCONFIG")
          << "deltaFunctionGun: ERROR: momentum magnitude can not be negative, got "<<pmag
          << "\n";
      }

      etot = sqrt(std::pow(pmag, 2) + std::pow(mass, 2));
    }
    else if(conf().kineticEnergy()) {
      double ek = *conf().kineticEnergy();
      if(ek < 0.) {
        throw cet::exception("BADCONFIG")
          << "deltaFunctionGun: ERROR: kinetic energy can not be negative, got "<<ek
          << "\n";
      }

      etot = ek + mass;
      pmag = sqrt(ek*(ek+2*mass));
    }
    else {
      throw cet::exception("BADCONFIG")
        << "deltaFunctionGun: ERROR: either momentum or kineticEnergy has to be specified"
        << "\n";
    }

    //----------------------------------------------------------------
    CLHEP::Hep3Vector dir(conf().direction());
    const double norm = dir.mag();
    if(norm > 0.) {
      dir = dir/norm;
    }
    else {
      throw cet::exception("BADCONFIG")
        << "deltaFunctionGun: ERROR: the direction vector must be non-zero"
        << "\n";
    }

    const CLHEP::Hep3Vector mom3 = pmag * dir;

    CLHEP::HepLorentzVector fourmom(etot, mom3);

    return GenParticle(pdg.id(),
                       GenId::particleGun,
                       conf().position(),
                       fourmom,
                       0.,
                       0.);
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::deltaFunctionGun)
