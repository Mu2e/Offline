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
      unsigned int nRejected_;
      unsigned int nTotal_;
      double showerEnergyCutoff_;
      int seed_;
      art::RandomNumberGenerator::base_engine_t&     engine_;
  };

  CryEventGenerator::CryEventGenerator(fhicl::ParameterSet const& pSet) :
    inputfile(pSet.get<std::string>("inputFile",
          "CRYEventGenerator/config/defaultCRYconfig.txt")),
    nRejected_(0), nTotal_(0), showerEnergyCutoff_(pSet.get<double>("showerEnergyCutoff", 1E6)),
    seed_( art::ServiceHandle<SeedService>()->getSeed() ),
    engine_(createEngine(seed_))
  {
    mf::LogInfo("CRYEventGenerator") << "Cutoff energy: " << showerEnergyCutoff_ << " MeV.";
    produces<GenParticleCollection>();
  }

  void CryEventGenerator::beginRun( art::Run &run){
    cryGen = std::make_unique<CosmicCRY>(run, SimpleConfig(inputfile), engine_);
  }

  void CryEventGenerator::produce(art::Event& evt) {
    std::unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
    bool belowECutoff = false;
    // keep clear + generate partiles untill everything is below cut off
    while (!belowECutoff){
      belowECutoff = true;
      nTotal_ ++;
      genParticles->clear();
      cryGen->generate(*genParticles);

      if (cryGen->getShowerSumEnergy() > showerEnergyCutoff_) {
        // mf::LogInfo("CRYEventGenerator") << "total E: " << cryGen->getShowerSumEnergy();
        // for (unsigned int i = 0; i < genParticles->size(); ++i) {
          // mf::LogInfo("CRYEventGenerator") << "nRejected_ " << nRejected_ << ", particle: " << genParticles->at(i) << genParticles->at(i).momentum().e();
        // }
        belowECutoff = false;
      }

      if (!belowECutoff) 
        nRejected_ ++;
    }

    evt.put(std::move(genParticles));
  }

  void CryEventGenerator::endRun(art::Run&){
    std::ostringstream oss;
    oss << "Total live time simulated in this run: " << cryGen->getLiveTime() << "\n";
    oss << "Total number of events: " << nTotal_ << "\n";
    oss << "Number of events rejected due to shower energy cutoff (" << showerEnergyCutoff_ << " MeV): " << nRejected_;
    mf::LogInfo("CRYEventGenerator") << oss.str();
  }

}


using mu2e::CryEventGenerator;
DEFINE_ART_MODULE(CryEventGenerator);
