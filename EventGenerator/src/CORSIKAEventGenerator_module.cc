// Mu2e includes.
#include "ConfigTools/inc/SimpleConfig.hh"
#include "ConfigTools/inc/requireUniqueKey.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"
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

#include "EventGenerator/inc/CosmicCORSIKA.hh"

    namespace mu2e {

  class CorsikaEventGenerator : public art::EDProducer {
    public:
      explicit CorsikaEventGenerator(fhicl::ParameterSet const& pSet);
      // Accept compiler written d'tor.  Modules are never moved or copied.
      virtual void produce (art::Event& e);
      virtual void beginRun(art::Run&   r);
      virtual void endRun(art::Run&   r);
    private:
      std::unique_ptr<CosmicCORSIKA> corsikaGen;
      std::string inputfile;
      int _seed;
      art::RandomNumberGenerator::base_engine_t &_engine;
  };

  CorsikaEventGenerator::CorsikaEventGenerator(fhicl::ParameterSet const &pSet) : EDProducer{pSet},
                                                                                  inputfile(pSet.get<std::string>("inputFile", "EventGenerator/config/defaultCORSIKAconfig.txt")),
                                                                                  _seed( art::ServiceHandle<SeedService>()->getSeed() ),
                                                                                  _engine(createEngine(_seed))
  {
    produces<GenParticleCollection>();
  }

  void CorsikaEventGenerator::beginRun( art::Run &run){
    corsikaGen = std::make_unique<CosmicCORSIKA>(run,
                                                 SimpleConfig(inputfile),
                                                 _engine);
  }

  void CorsikaEventGenerator::produce(art::Event &evt)
  {
    std::unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
    genParticles->clear();
    corsikaGen->generate(*genParticles);
    evt.put(std::move(genParticles));
  }

  void CorsikaEventGenerator::endRun(art::Run &)
  {
    std::ostringstream oss;
    mf::LogInfo("CORSIKAEventGenerator") << "Total primaries: " << corsikaGen->getNumShowers(); // << std::endl;
    mf::LogInfo("CORSIKAEventGenerator") << "Total live-time: " << corsikaGen->getLiveTime(); // << std::endl;

    mf::LogInfo("CORSIKAEventGenerator") << oss.str();
  }
}


using mu2e::CorsikaEventGenerator;
DEFINE_ART_MODULE(CorsikaEventGenerator);
