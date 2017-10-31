// Mu2e includes.
#include "ConfigTools/inc/SimpleConfig.hh"
#include "ConfigTools/inc/requireUniqueKey.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
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

// Other external includes.
#include <boost/shared_ptr.hpp>

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CRYEventGenerator/inc/CosmicCRY.hh"

namespace mu2e {

  class CryEventGenerator : public art::EDProducer {
    public:
      explicit CryEventGenerator(fhicl::ParameterSet const& pSet);
      // Accept compiler written d'tor.  Modules are never moved or copied.
      virtual void produce (art::Event& e);
      virtual void beginRun(art::Run&   r);
    private:
      CosmicCRY *cryGen;
      std::string inputfile;
  };

  CryEventGenerator::CryEventGenerator(fhicl::ParameterSet const& pSet) :
    inputfile(pSet.get<std::string>("inputFile",
          "CRYEventGenerator/config/defaultCRYconfig.txt"))
  {
    // createEngine(art::ServiceHandle<SeedService>()->getSeed());
    
    // use G4Engine to correctly get seed from fcl file
    createEngine(art::ServiceHandle<SeedService>()->getSeed(), "G4Engine");
    produces<GenParticleCollection>();
  }

  void CryEventGenerator::beginRun( art::Run &run){
    cryGen = new CosmicCRY(run, SimpleConfig(inputfile));
  }

  void CryEventGenerator::produce(art::Event& evt) {
    std::unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
    cryGen->generate(*genParticles);
    evt.put(std::move(genParticles));
  }


}


using mu2e::CryEventGenerator;
DEFINE_ART_MODULE(CryEventGenerator);
