// Sophie Middleton, 2022
// C++ includes.
#include <iostream>
#include <algorithm>

// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "Offline/ConditionsService/inc/AcceleratorParams.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataTable.hh"
//#include "Offline/EventGenerator/inc/CaloCalibGun.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"

// Other external includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

//ROOT Includes
#include "TMath.h"
using namespace std;
using namespace TMath;

// Generates particles in a flat shaped spectrum, size determined by FCL params these will be attached to a mu- in
// the input SimParticleCollection.
// This module throws an exception if no suitable muon is found.
//
// S Middleton, 2021

#include <iostream>
#include <string>
#include <cmath>
#include <memory>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"

//For primary:
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"

#include <TTree.h>
namespace mu2e {

  //================================================================
  class NewCaloCalibGun : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<double> PhotonEnergy{Name("PhotonEnergy"),6.13};//MeV
      fhicl::Atom<double> Mean{Name("Mean"),1.},
      fhicl::Atom<double> CosMin{Name("cosmin"),-1.},
      fhicl::Atom<double> CosMin{Name("cosmin"),1.},
      fhicl::Atom<double> PhiMin{Name("phimin"),0.},
      fhicl::Atom<double> PhiMin{Name("phimax"), CLHEP::twopi },
      fhicl::Atom<bool> doHistograms{Name("doHistograms"), true },
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit NewCaloCalibGun(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:

    double PhotonEnergy_;
    
  };

  //================================================================
  NewCaloCalibGun::NewCaloCalibGun(const Parameters& conf)
    : EDProducer{conf}
    , PhotonEnergy_{conf().PhotonEnergy()}
    , randFlat_{engine},
    , randPoissonQ_{engine, std::abs(_mean)},
    , randomUnitSphere_{engine, _cosmin, _cosmax, 0, CLHEP::twopi},
//  _detSys(),
  {
    produces<mu2e::StageParticleCollection>();
    produces<mu2e::PrimaryParticle>();
  }

  //================================================================
  void NewCaloCalibGun::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};
    //Call the gun here
    output->emplace_back(mustop,
                         ProcessCode::CaloCalib,
                         PDGCode::gamma,
                         mustop->endPosition(),
                         CLHEP::HepLorentzVector{randomUnitSphere_.fire(randomMom), randomE},
                         time
                         );
    
    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::NewCaloCalibGun);
