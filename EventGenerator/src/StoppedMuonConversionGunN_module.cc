// Modification of StopppedMuonConversionGun to act as a source of particles
// for the Tom LeCompte style G4MT mode.
//
// Rob Kutschke, 2017

#include <string>
#include <cmath>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollections.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"

namespace mu2e {

  //================================================================
  class StoppedMuonConversionGunN : public art::EDProducer {
    double conversionEnergy_;
    double conversionMomentum_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    double           czmin_;
    double           czmax_;
    double           phimin_;
    double           phimax_;
    RandomUnitSphere randomUnitSphere_;

    RootTreeSampler<IO::StoppedParticleF> stops_;

    size_t stashSize_;
    size_t finger_;
    GenParticleCollections stash_;

    static double electronMass() {
      return GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::e_minus).ref().mass().value();
    }

  public:
    explicit StoppedMuonConversionGunN(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
    void resetStash();
    void fillStash();
  };

  //================================================================
  StoppedMuonConversionGunN::StoppedMuonConversionGunN(const fhicl::ParameterSet& pset)
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
    , stops_           (eng_, pset.get<fhicl::ParameterSet>("muonStops"))
    , stashSize_       (pset.get<size_t>("stashSize"))
    , finger_          (stashSize_)
    , stash_           (stashSize_)
  {
    produces<mu2e::GenParticleCollection>();
    produces<mu2e::GenParticleCollections>();
  }

  //================================================================
  void StoppedMuonConversionGunN::produce(art::Event& event) {

    // On event 0, N, 2N, ... fill the stash with N events
    // and place a copy of the stash in the event.
    if ( finger_ == stash_.size() ){
      resetStash();
      fillStash();
      finger_=0;
      auto allEvents = std::make_unique<GenParticleCollections>(stash_);
      event.put(std::move(allEvents));
    } else{
      auto allEvents = std::make_unique<GenParticleCollections>();
      event.put(std::move(allEvents));
    }

    // On every event, copy the next GenParticleCollection out of the stash and
    // put it into the event.
    auto gens = std::make_unique<GenParticleCollection>(stash_[finger_++]);
    event.put(std::move(gens));
  }

  void StoppedMuonConversionGunN::resetStash() {
    for ( size_t i=0; i<stashSize_; ++i){
      stash_[i].clear();
    }
  }

  void StoppedMuonConversionGunN::fillStash() {

    for ( size_t i=0; i<stashSize_; ++i){

      const auto& stop = stops_.fire();
      const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

      const CLHEP::Hep3Vector p3 = randomUnitSphere_.fire(conversionMomentum_);
      const CLHEP::HepLorentzVector p4(p3, conversionEnergy_);

      stash_[i].emplace_back(PDGCode::e_minus, GenId::conversionGun, pos, p4, stop.t);
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedMuonConversionGunN);
