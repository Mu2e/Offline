// Lisa Goodenough, 2018
/*
 
Modification of StoppedParticleReactionGun for running in MT mode.
It produces a collection of GenParticleCollections which
are used in the stashes in MT mode of Mu2eG4_module.cc.
*/


#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollections.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/EjectedProtonSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"

namespace mu2e {

  //================================================================
  class StoppedParticleReactionGunN : public art::EDProducer {
    fhicl::ParameterSet psphys_;

    PDGCode::type pdgId_;
    double mass_;

    enum SpectrumVar { TOTAL_ENERGY, KINETIC_ENERY, MOMENTUM };
    SpectrumVar spectrumVariable_;
    static SpectrumVar parseSpectrumVar(const std::string& name);

    BinnedSpectrum spectrum_;
    
    GenId genId_;
      
    int verbosityLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandGeneral randSpectrum_;
    RandomUnitSphere randomUnitSphere_;

    RootTreeSampler<IO::StoppedParticleF> stops_;
      
    size_t stashSize_;
    size_t finger_;
    GenParticleCollections stash_;

    double generateEnergy();

  public:
    explicit StoppedParticleReactionGunN(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
    void resetStash();
    void fillStash();
  };

  //================================================================
  StoppedParticleReactionGunN::StoppedParticleReactionGunN(const fhicl::ParameterSet& pset)
    : art::EDProducer{pset}
    , psphys_(pset.get<fhicl::ParameterSet>("physics"))
    , pdgId_(PDGCode::type(psphys_.get<int>("pdgId")))
    , mass_(GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId_).ref().mass().value())
    , spectrumVariable_(parseSpectrumVar(psphys_.get<std::string>("spectrumVariable")))
    , spectrum_(BinnedSpectrum(psphys_))
    , genId_(GenId::findByName(psphys_.get<std::string>("genId", "StoppedParticleReactionGun")))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randSpectrum_(eng_, spectrum_.getPDF(), spectrum_.getNbins())
    , randomUnitSphere_(eng_)
    , stops_(eng_, pset.get<fhicl::ParameterSet>("muonStops"))
    , stashSize_       (pset.get<size_t>("stashSize"))
    , finger_          (stashSize_)
    , stash_           (stashSize_)
  {
    produces<mu2e::GenParticleCollection>();
    produces<mu2e::GenParticleCollections>();
      
    if(genId_ == GenId::enum_type::unknown) {
          throw cet::exception("BADCONFIG")<<"StoppedParticleReactionGun: unknown genId "
          <<psphys_.get<std::string>("genId", "StoppedParticleReactionGun")
          <<"\n";
    }

    if(verbosityLevel_ > 0) {
      std::cout<<"StoppedParticleReactionGunN: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;

      std::cout<<"StoppedParticleReactionGun: using GenId = "
               <<genId_
               <<std::endl;
        
      std::cout<<"StoppedParticleReactionGunN: producing particle "
               <<pdgId_
               <<", mass = "<<mass_
               <<std::endl;

      std::cout <<"StoppedParticleReactionGunN: spectrum shape = "
	  <<psphys_.get<std::string>("spectrumShape") << std::endl;
      if (psphys_.get<std::string>("spectrumShape")  == "tabulated")
	  std::cout << " Spectrum file = "
	  << psphys_.get<std::string>("spectrumFileName")
	  << std::endl;
    }
    if(verbosityLevel_ > 1){
      std::cout <<"StoppedParticleReactionGunN: spectrum: " << std::endl;
      spectrum_.print();
    }
  }

  //================================================================
  StoppedParticleReactionGunN::SpectrumVar
  StoppedParticleReactionGunN::parseSpectrumVar(const std::string& name) {
    if(name == "totalEnergy")  return TOTAL_ENERGY;
    if(name == "kineticEnergy")  return KINETIC_ENERY;
    if(name == "momentum")  return MOMENTUM;
    throw cet::exception("BADCONFIG")<<"StoppedParticleReactionGunN: unknown spectrum variable "<<name<<"\n";
  }


  //================================================================
  void StoppedParticleReactionGunN::produce(art::Event& event) {

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
    
    //================================================================
  void StoppedParticleReactionGunN::resetStash() {
    for ( size_t i=0; i<stashSize_; ++i){
        stash_[i].clear();
    }
  }
    
  //================================================================
  void StoppedParticleReactionGunN::fillStash() {
        
    for ( size_t i=0; i<stashSize_; ++i){
        const auto& stop = stops_.fire();
        const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
        
        const double energy = generateEnergy();
        const double p = energy * sqrt(1 - std::pow(mass_/energy,2));
        
        CLHEP::Hep3Vector p3 = randomUnitSphere_.fire(p);
        CLHEP::HepLorentzVector fourmom(p3, energy);
        
        stash_[i].emplace_back(pdgId_,
                               genId_,
                               pos,
                               fourmom,
                               stop.t);
        
    }
  }

  //================================================================
  double StoppedParticleReactionGunN::generateEnergy() {
    double res = spectrum_.sample(randSpectrum_.fire());
    if(res < 0.0)
        throw cet::exception("BADE")<<"StoppedParticleReactionGun: negative energy "<< res <<"\n";
    switch(spectrumVariable_) {
    case TOTAL_ENERGY: break;
    case KINETIC_ENERY: res += mass_; break;
    case MOMENTUM     : res = sqrt(res*res+mass_*mass_); break;
    }
    return res;
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedParticleReactionGunN);
