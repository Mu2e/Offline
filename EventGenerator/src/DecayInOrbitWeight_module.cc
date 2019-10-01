// Kyle Knoepfel, 2014

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
#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Mu2eUtilities/inc/ShankerWatanabeSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"

namespace mu2e {

  //================================================================
  class DecayInOrbitWeight : public art::EDProducer {

    class SpectrumType {
    public:
      enum enum_type { unknown, Pol5, Pol58 };
      static std::string const& typeName() {
        static std::string type("SpectrumType"); return type;
      }
      static std::map<enum_type,std::string> const& names() {
        static std::map<enum_type,std::string> nam;
        
        if ( nam.empty() ) {
          nam[unknown] = "unknown";
          nam[Pol5]    = "pol5";
          nam[Pol58]   = "pol58";
        }

        return nam;
      }
    };

    typedef EnumToStringSparse<SpectrumType> SpectrumChoice;

    art::InputTag input_;
    SpectrumChoice weightingScheme_;

    int verbosityLevel_;
  public:
    explicit DecayInOrbitWeight(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
  };

  //================================================================
  DecayInOrbitWeight::DecayInOrbitWeight(const fhicl::ParameterSet& pset)
    : art::EDProducer{pset}
    , input_(pset.get<std::string>("inputModule") )
    , weightingScheme_(pset.get<std::string>("weightingScheme","pol58") )
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0 ) )
  {
    produces<mu2e::EventWeight>();

    if(verbosityLevel_ > 0) {
      std::cout<< "DecayInOrbitWeight: using the "
               << weightingScheme_
               << " weighting scheme"
               << std::endl;

      std::cout<< "DecayInOrbitWeight: "
               << std::endl << std::endl
               << " NOTE: " 
               << std::endl << std::endl
               << " The event weights added in this module should only be used when " << std::endl
               << " the corresponding spectrum in DecayInOrbitGun is 'flat'. " 
               << std::endl;
    }
  }

  //================================================================
  void DecayInOrbitWeight::produce(art::Event& event) {

    auto genColl = event.getValidHandle<GenParticleCollection>( input_ );

    double weight(-1.);
    for ( const auto& i: *genColl ) {
      if (i.generatorId() == GenId::dioTail) {
	const double energy = i.momentum().e();
      
	if      ( weightingScheme_ == SpectrumChoice::Pol5  ) weight = SimpleSpectrum::getPol5 ( energy ); 
	else if ( weightingScheme_ == SpectrumChoice::Pol58 ) weight = SimpleSpectrum::getPol58( energy ); 
	else {
	  throw cet::exception("MODEL")
	    << "Wrong or not allowed DIO energy spectrum";
	}
      }
    }

    std::unique_ptr<mu2e::EventWeight> evtwt ( new EventWeight(weight) );
    event.put(std::move(evtwt));
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::DecayInOrbitWeight);
