// Andrei Gaponenko, 2014

#include <string>
#include <cmath>

#include "cetlib/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "SeedService/inc/SeedService.hh"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/PhysicsParams.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
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
  }

  //================================================================
  class StoppedMuonConversionGun : public art::EDProducer {
    double conversionEnergy_;
    double conversionMomentum_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    double           czmin_;
    double           czmax_;
    double           phimin_;
    double           phimax_;
    RandomUnitSphere randomUnitSphere_;

    RootTreeSampler<InputStop> stops_;

    static double electronMass() {
      return GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::e_minus).ref().mass().value();
    }

  public:
    explicit StoppedMuonConversionGun(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
  };

  //================================================================
  StoppedMuonConversionGun::StoppedMuonConversionGun(const fhicl::ParameterSet& pset)
    : conversionEnergy_(GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy())
    , conversionMomentum_(conversionEnergy_*
                          sqrt(1 - std::pow(electronMass()/conversionEnergy_,2))
                          )
    , eng_             (createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , czmin_           (pset.get<double>("czmin" , -1.0))
    , czmax_           (pset.get<double>("czmax" ,  1.0))
    , phimin_          (pset.get<double>("phimin",  0. ))
    , phimax_          (pset.get<double>("phimax", CLHEP::twopi ))
    , randomUnitSphere_(eng_,czmin_,czmax_,phimin_,phimax_)
    , stops_           (eng_, pset.get<fhicl::ParameterSet>("muonStops"), sizeof(InputStop)/sizeof(float))
  {
    produces<mu2e::GenParticleCollection>();
  }

  //================================================================
  void StoppedMuonConversionGun::produce(art::Event& event) {
    const InputStop& stop = stops_.fire();
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    const CLHEP::Hep3Vector p3 = randomUnitSphere_.fire(conversionMomentum_);
    const CLHEP::HepLorentzVector p4(p3, conversionEnergy_);

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    output->emplace_back(PDGCode::e_minus, GenId::conversionGun, pos, p4, stop.t);
    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedMuonConversionGun);
