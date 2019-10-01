// Andy Edmonds, 2018

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
#include "MCDataProducts/inc/EventWeight.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"

namespace mu2e {

  //================================================================
  class PBIWeight : public art::EDProducer {

    art::InputTag _PBITag;
    art::InputTag _meanPBITag;
    double _meanPBI;
  public:
    explicit PBIWeight(const fhicl::ParameterSet& pset);
    virtual void beginSubRun(art::SubRun& subrun ) override;
    virtual void produce(art::Event& event);
  };

  //================================================================
  PBIWeight::PBIWeight(const fhicl::ParameterSet& pset) : 
    art::EDProducer{pset},
    _PBITag(pset.get<art::InputTag>("PBITag")),
    _meanPBITag(pset.get<art::InputTag>("meanPBITag")),
    _meanPBI(0.0)
	    
  {
    produces<mu2e::EventWeight>();
  }

  void PBIWeight::beginSubRun(art::SubRun& subrun) {
    art::Handle<ProtonBunchIntensity> PBIHandle;
    subrun.getByLabel(_meanPBITag, PBIHandle);
    if(PBIHandle.isValid()) {
      _meanPBI = PBIHandle->intensity();
    }
  }

  //================================================================
  void PBIWeight::produce(art::Event& event) {

    art::Handle<ProtonBunchIntensity> PBIHandle;
    event.getByLabel(_PBITag, PBIHandle);
    double weight = 1.0; // by default weight will be 1 (e.g. for non-background mixed jobs)
    if (PBIHandle.isValid()) {
      double PBI = PBIHandle->intensity();
      weight = PBI / _meanPBI;
    }
    std::unique_ptr<mu2e::EventWeight> evtwt ( new EventWeight(weight) );
    event.put(std::move(evtwt));
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::PBIWeight);
