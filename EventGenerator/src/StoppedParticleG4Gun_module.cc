// This generator places particle of a given pdg ID at points sampled from
// an input ROOT tree.  The particles are produced at rest.
//
// One can not replace this generator with a special case of
// StoppedParticleReactionGun, because the timing of the produced
// particles needs different handling downstream: a random delay due
// to muon life time should be added to StoppedParticleReactionGun
// products (for daughters of stopped muons), but not to G4-produced
// daughters of muons that were placed by this generator.
//
// Andrei Gaponenko, 2014

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
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/Mu2eUtilities/inc/RootTreeSampler.hh"
#include "Offline/GeneralUtilities/inc/RSNTIO.hh"

namespace mu2e {

  //================================================================
  class StoppedParticleG4Gun : public art::EDProducer {
    PDGCode::type pdgId_;
    double mass_;
    int verbosityLevel_;

    typedef RootTreeSampler<IO::StoppedParticleF> RTS;
    RTS stops_;

  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<int> pdgId { Name("pdgId"), Comment("PDG ID of the particle to place.")};

      fhicl::Atom<int> verbosityLevel{
        Name("verbosityLevel"),
          Comment("A positive value generates more printouts than the default."),
          0
          };

      fhicl::Table<RTS::Config> muonStops {
        Name("muonStops"),
          Comment("RootTreeSampler settings")
          };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit StoppedParticleG4Gun(const Parameters& conf);

    virtual void produce(art::Event& event);
  };

  //================================================================
  StoppedParticleG4Gun::StoppedParticleG4Gun(const Parameters& conf)
    : EDProducer{conf}
    , pdgId_(PDGCode::type(conf().pdgId()))
    , mass_(GlobalConstantsHandle<ParticleDataList>()->particle(pdgId_).mass())
    , verbosityLevel_(conf().verbosityLevel())
    , stops_(createEngine(art::ServiceHandle<SeedService>()->getSeed()),
             conf().muonStops())
  {
    produces<mu2e::GenParticleCollection>();

    if(verbosityLevel_ > 0) {
      std::cout<<"StoppedParticleG4Gun: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;

      std::cout<<"StoppedParticleG4Gun: producing particle "
               <<pdgId_
               <<", mass = "<<mass_
               <<std::endl;
    }
  }

  //================================================================
  void StoppedParticleG4Gun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    const auto& stop = stops_.fire();
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    CLHEP::HepLorentzVector fourmom(CLHEP::Hep3Vector(), mass_);

    output->emplace_back(pdgId_,
                         GenId::StoppedParticleG4Gun,
                         pos,
                         fourmom,
                         stop.t);

    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedParticleG4Gun)
