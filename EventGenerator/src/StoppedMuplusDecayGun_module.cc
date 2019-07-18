// Andrei Gaponenko, 2014
// adapted from StoppedMuoneDecayGun by D. Brown, Aug. 2015
// Decays stopped muons according to pure Michel spectrum,
// without an nuclear binding (ie only correct for mu+ decay)

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

#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"
#include "CLHEP/Random/RandFlat.h"

namespace mu2e {

  //================================================================
  class StoppedMuplusDecayGun : public art::EDProducer {

    art::RandomNumberGenerator::base_engine_t& eng_;
    double           czmin_;
    double           czmax_;
    double           phimin_;
    double           phimax_;
    double           emin_; // minimum electron energy
    RandomUnitSphere randomUnitSphere_;
    RootTreeSampler<IO::StoppedParticleF> stops_;
    CLHEP::RandFlat flat_;

    static double electronMass2() {
      static double emass = GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::e_minus).ref().mass().value();
      return emass*emass;
    }
    static double muMass() {
      static double mumass = GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::mu_minus).ref().mass().value();
      return mumass;
    }

  public:
    explicit StoppedMuplusDecayGun(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
  };

  //================================================================
  StoppedMuplusDecayGun::StoppedMuplusDecayGun(const fhicl::ParameterSet& pset) :
    EDProducer{pset},
    eng_             (createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , czmin_           (pset.get<double>("czmin" , -1.0))
    , czmax_           (pset.get<double>("czmax" ,  1.0))
    , phimin_          (pset.get<double>("phimin",  0. ))
    , phimax_          (pset.get<double>("phimax", CLHEP::twopi ))
    , emin_            (pset.get<double>("emin", 10.0))
    , randomUnitSphere_(eng_,czmin_,czmax_,phimin_,phimax_)
    , stops_           (eng_, pset.get<fhicl::ParameterSet>("muonStops"))
    , flat_        (eng_)
  {
    produces<mu2e::GenParticleCollection>();
  }

  //================================================================
  void StoppedMuplusDecayGun::produce(art::Event& event) {
    const auto& stop = stops_.fire();
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

// generate the electron momentum.  Just use throw/miss since the spectrum is simple
    bool pass(false);
    double xval(emin_);
    while(!pass){
    // xval = ratio of energy to half the muon mass
      xval = flat_.fire(2.0*emin_/muMass(),1.0);
      double spect = xval*xval*(3.0-2.0*xval);
      pass = flat_.fire(0.0,1.0) < spect;
    }
    double eele = xval*0.5*muMass();
    double pele = sqrt(eele*eele - electronMass2());
    const CLHEP::Hep3Vector p3 = randomUnitSphere_.fire(pele);
    const CLHEP::HepLorentzVector p4(p3, eele);

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    output->emplace_back(PDGCode::e_plus, GenId::muplusDecayGun, pos, p4, stop.t);
    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedMuplusDecayGun);
