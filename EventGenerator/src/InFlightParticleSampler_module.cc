// This generator samples particle position, momentum, and pdg ID from
// a ROOT tree.  Used to re-sample particles intercepted by the
// "killer planes" in a previous simulation stage.
//
// Andrei Gaponenko, 2015

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"

namespace mu2e {

  //================================================================
  class InFlightParticleSampler : public art::EDProducer {
    typedef RootTreeSampler<std::vector<IO::InFlightParticleD>, IO::InFlightParticleD > RTS;
    RTS particles_;

    const mu2e::ParticleDataTable *pdt_;
    int verbosityLevel_;

  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<int> verbosityLevel{
        Name("verbosityLevel"),
          Comment("A positive value generates more printouts than the default."),
          0
          };

      fhicl::Table<RTS::Config> particles {
        Name("particles"),
          Comment("RootTreeSampler settings")
          };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit InFlightParticleSampler(const Parameters& conf);

    virtual void produce(art::Event& event);
  };
  //================================================================
  InFlightParticleSampler::InFlightParticleSampler(const Parameters& conf)
    : EDProducer{conf}
    , particles_(createEngine(art::ServiceHandle<SeedService>()->getSeed()),
                 conf().particles())

    , pdt_(&*GlobalConstantsHandle<ParticleDataTable>())

    , verbosityLevel_(conf().verbosityLevel())
  {
    produces<mu2e::GenParticleCollection>();

    if(verbosityLevel_ > 0) {
      std::cout<<"InFlightParticleSampler: using = "
               <<particles_.numRecords()
               <<" particles"
               <<std::endl;
    }
  }

  //================================================================
  void InFlightParticleSampler::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    const auto& list = particles_.fire();

    for(const auto& part: list) {
      const CLHEP::Hep3Vector p3(part.px, part.py, part.pz);
      const double mass = pdt_->particle(part.pdgId).ref().mass().value();
      const double energy = std::sqrt(std::pow(mass,2) + p3.mag2());
      CLHEP::HepLorentzVector fourmom(p3, energy);

      const CLHEP::Hep3Vector pos(part.x, part.y, part.z);

      output->emplace_back(PDGCode::type(part.pdgId),
                           GenId::InFlightParticleSampler,
                           pos,
                           fourmom,
                           part.time);
    }

    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::InFlightParticleSampler);
