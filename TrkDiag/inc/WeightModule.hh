//Hasung Song (2018)

// C++ includes
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

// cetlib includes
#include "cetlib_except/exception.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

// Mu2e includes
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"

namespace mu2e {

  template <class Phys> 

  class WeightModule : public art::EDProducer {
    Phys _myPhys;
    public:
    explicit WeightModule(const fhicl::ParameterSet& pset):
      art::EDProducer{pset},
      _myPhys(pset)
    {
      produces<mu2e::EventWeight>();
    } ;
    virtual void produce(art::Event& event) {
      double weight (-1.);
      weight = _myPhys.weight(event);
      std::unique_ptr<mu2e::EventWeight> evtwt ( new EventWeight(weight) );
      event.put(std::move(evtwt));
    };

  };

} //namespace mu2e
