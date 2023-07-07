//
// A generator to shoot electrons right in the face of the caloriemter disk to simulate the calo test beam. Fire!
//

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// Mu2e includes
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

#include <cmath>

namespace mu2e {

  //================================================================
  class CaloTBGun : public art::EDProducer
  {
     public:
       explicit     CaloTBGun  (const fhicl::ParameterSet& pset);
       virtual void produce (art::Event& event);

     private:
       int                 verbosityLevel_;
       double              energy_;
       double              angle_;
       double              time_;
       bool                frontVD_;
  };

  //================================================================
  CaloTBGun::CaloTBGun(const fhicl::ParameterSet& pset)
    : EDProducer{pset}
    , verbosityLevel_            (pset.get<int>   ("verbosityLevel", 0    ))
    , energy_                    (pset.get<double>("energy",         105.0))
    , angle_                     (pset.get<double>("angle",          0.0  ))
    , time_                      (pset.get<double>("time",           700.0))
    , frontVD_                   (pset.get<bool>  ("frontVD",        true ))

  {
      produces<mu2e::GenParticleCollection>();
      if (verbosityLevel_ > 0) std::cout<<"CaloTB gun: shoot! " << std::endl;
  }

  //================================================================
  void CaloTBGun::produce(art::Event& event)
  {
       const Calorimeter& cal = *(GeomHandle<Calorimeter>());

       double                  angleRad = angle_*M_PI/180;
       double                  Zbuffer  = cal.disk(0).crystal(400).position().z()-cal.disk(0).geomInfo().origin().z()+cal.disk(0).geomInfo().size().z()/2.0+1;
       double                  dz       = (frontVD_)? Zbuffer : 0;
       CLHEP::Hep3Vector       pos      = cal.disk(0).crystal(400).position() - CLHEP::Hep3Vector(tan(angleRad)*dz,0,dz);
       CLHEP::HepLorentzVector mom(energy_*sin(angleRad),0.0,energy_*cos(angleRad),energy_);

       std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
       output->emplace_back(PDGCode::e_minus,GenId::CeEndpoint,pos,mom,time_);
       event.put(std::move(output));
  }

}

DEFINE_ART_MODULE(mu2e::CaloTBGun)
