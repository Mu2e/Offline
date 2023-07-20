// Mu2e includes.
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/ConfigTools/inc/requireUniqueKey.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/G4BeamlineInfo.hh"
#include "Offline/MCDataProducts/inc/CosmicLivetime.hh"

// Particular generators that this code knows about.
#include "Offline/SeedService/inc/SeedService.hh"

// Includes from art and its toolchain.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>

#include "Offline/EventGenerator/inc/CosmicCRY.hh"

namespace mu2e {

  class CryEventGenerator : public art::EDProducer {
    public:
      explicit CryEventGenerator(fhicl::ParameterSet const& pSet);
      // Accept compiler written d'tor.  Modules are never moved or copied.
      virtual void produce (art::Event& e);
      virtual void beginRun(art::Run&   r);
      virtual void endRun(art::Run&   r);
      virtual void endSubRun(art::SubRun &sr);
    private:
      std::unique_ptr<CosmicCRY> cryGen;
      std::string inputfile;
      art::RandomNumberGenerator::base_engine_t&     engine_;
  };

  CryEventGenerator::CryEventGenerator(fhicl::ParameterSet const& pSet) :
    EDProducer{pSet},
    inputfile(pSet.get<std::string>("inputFile",
          "CRYEventGenerator/config/defaultCRYconfig.txt")),
    engine_( createEngine( art::ServiceHandle<SeedService>()->getSeed()) )
  {
    produces<GenParticleCollection>();
    produces<CosmicLivetime,art::InSubRun>();
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

  void CryEventGenerator::endSubRun(art::SubRun &subrun)
  {
    // All inputs (except getLiveTime) in CosmicLivetime are not used in the livetime calculation
    // Livetime is calculated internally by cry: cryGen->getLiveTime()
    std::unique_ptr<CosmicLivetime> livetime(new CosmicLivetime(1,
                                                                cryGen->getSubboxLength()*cryGen->getSubboxLength(),
                                                                cryGen->getMinShowerEn()/1000.,
                                                                cryGen->getMaxShowerEn()/1000.,
                                                                1.8e4, // http://pdg.lbl.gov/2018/reviews/rpp2018-rev-cosmic-rays.pdf eq. 29.2
                                                                cryGen->getLiveTime()  ));
    std::cout << *livetime << std::endl;
    subrun.put(std::move(livetime), art::fullSubRun());
  }

}


using mu2e::CryEventGenerator;
DEFINE_ART_MODULE(CryEventGenerator)
