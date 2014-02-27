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

#include "cetlib/exception.h"

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
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/PhysicsParams.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"

namespace mu2e {

  //================================================================
  namespace {
    struct InputStop {
      float x;
      float y;
      float z;
      float t;

      InputStop() : x(), y(), z(), t() {}
    };
  } // namespace {}

  //================================================================
  class StoppedParticleG4Gun : public art::EDProducer {
    PDGCode::type pdgId_;
    double mass_;
    int verbosityLevel_;

    RootTreeSampler<InputStop> stops_;

  public:
    explicit StoppedParticleG4Gun(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
  };

  //================================================================
  StoppedParticleG4Gun::StoppedParticleG4Gun(const fhicl::ParameterSet& pset)
    : pdgId_(PDGCode::type(pset.get<int>("pdgId")))
    , mass_(GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId_).ref().mass().value())
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , stops_(createEngine(art::ServiceHandle<SeedService>()->getSeed()),
             pset.get<fhicl::ParameterSet>("muonStops"), sizeof(InputStop)/sizeof(float))
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

    const InputStop& stop = stops_.fire();
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

DEFINE_ART_MODULE(mu2e::StoppedParticleG4Gun);
