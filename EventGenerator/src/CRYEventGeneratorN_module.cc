// Mu2e includes.
#include "ConfigTools/inc/SimpleConfig.hh"
#include "ConfigTools/inc/requireUniqueKey.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollections.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"

// Particular generators that this code knows about.
#include "SeedService/inc/SeedService.hh"

// Includes from art and its toolchain.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>

#include "EventGenerator/inc/CosmicCRY.hh"

namespace mu2e {

  class CryEventGeneratorN : public art::EDProducer {
    public:
      explicit CryEventGeneratorN(fhicl::ParameterSet const& pSet);
      // Accept compiler written d'tor.  Modules are never moved or copied.
      virtual void produce (art::Event& e);
      virtual void beginRun(art::Run&   r);
      virtual void endRun(art::Run&   r);
      void resetStash();
      void fillStash();

    private:
      std::unique_ptr<CosmicCRY> cryGen;
      std::string inputfile;
      int seed_;
      art::RandomNumberGenerator::base_engine_t&     engine_;
      size_t stashSize_;
      size_t finger_;
      GenParticleCollections stash_;

  };

  CryEventGeneratorN::CryEventGeneratorN(fhicl::ParameterSet const& pSet) :
    art::EDProducer{pSet},
    inputfile(pSet.get<std::string>("inputFile",
          "CRYEventGenerator/config/defaultCRYconfig.txt")),
    seed_( art::ServiceHandle<SeedService>()->getSeed() ),
    engine_(createEngine(seed_)),
    stashSize_(             pSet.get<size_t>        ("stashSize",             1)),
    finger_(stashSize_),
    stash_(stashSize_)    

  {
    produces<GenParticleCollection>();
    produces<GenParticleCollections>();
  }

  void CryEventGeneratorN::beginRun( art::Run &run){
    cryGen = std::make_unique<CosmicCRY>(run, SimpleConfig(inputfile), engine_);
  }

  void CryEventGeneratorN::produce(art::Event& evt) {
    //*****************************************************************
    // On event 0, N, 2N, ... fill the stash with N events
    // and place a copy of the stash in the event.
    if ( finger_ == stash_.size() ){
      resetStash();
      fillStash();
      finger_=0;
      auto allEvents = std::make_unique<GenParticleCollections>(stash_);
      evt.put(std::move(allEvents));
    } else{
      auto allEvents = std::make_unique<GenParticleCollections>();
      evt.put(std::move(allEvents));
    }
    // On every event, copy the next GenParticleCollection out of the stash and
    // put it into the event.
    auto gens = std::make_unique<GenParticleCollection>(stash_[finger_++]);
    evt.put(std::move(gens));
  }

  void CryEventGeneratorN::resetStash() {
    for ( size_t i=0; i<stashSize_; ++i){
      stash_[i].clear();
    }
  }
    
  void CryEventGeneratorN::fillStash() {    
    for ( size_t i=0; i<stashSize_; ++i){
      cryGen->generate(stash_[i]);
    }
  }

  void CryEventGeneratorN::endRun(art::Run&){
    std::ostringstream oss;
    oss << "Total live time simulated: " << cryGen->getLiveTime() << "\n";
    oss << "Number of events simulated: " << cryGen->getNumEvents() << "\n";
    mf::LogInfo("CRYEventGenerator") << oss.str();
  }

}


using mu2e::CryEventGeneratorN;
DEFINE_ART_MODULE(CryEventGeneratorN);
