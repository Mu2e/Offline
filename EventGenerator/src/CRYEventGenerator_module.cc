// Mu2e includes.
#include "ConfigTools/inc/SimpleConfig.hh"
#include "ConfigTools/inc/requireUniqueKey.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

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

  class CryEventGenerator : public art::EDProducer {
    public:
      explicit CryEventGenerator(fhicl::ParameterSet const& pSet);
      // Accept compiler written d'tor.  Modules are never moved or copied.
      virtual void produce (art::Event& e);
      virtual void beginRun(art::Run&   r);
      virtual void endRun(art::Run&   r);
    private:
      std::unique_ptr<CosmicCRY> cryGen;
      std::string inputfile;
      int seed_;
      art::RandomNumberGenerator::base_engine_t&     engine_;
  };

  CryEventGenerator::CryEventGenerator(fhicl::ParameterSet const& pSet) :
    EDProducer{pSet},
    inputfile(pSet.get<std::string>("inputFile",
          "CRYEventGenerator/config/defaultCRYconfig.txt")),
    seed_( art::ServiceHandle<SeedService>()->getSeed() ),
    engine_(createEngine(seed_))
  {
    produces<GenParticleCollection>();
  }

  void CryEventGenerator::beginRun( art::Run &run){
    cryGen = std::make_unique<CosmicCRY>(run, SimpleConfig(inputfile), engine_);
  }

  void CryEventGenerator::produce(art::Event& evt) {
    std::unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
    genParticles->clear();
    cryGen->generate(*genParticles);
    evt.put(std::move(genParticles));
  }

  void CryEventGenerator::endRun(art::Run&){
    std::ostringstream oss;
    oss << "Total live time simulated: " << cryGen->getLiveTime() << "\n";
    oss << "Number of events simulated: " << cryGen->getNumEvents() << "\n";
    mf::LogInfo("CRYEventGenerator") << oss.str();
  }

}


using mu2e::CryEventGenerator;
DEFINE_ART_MODULE(CryEventGenerator);
