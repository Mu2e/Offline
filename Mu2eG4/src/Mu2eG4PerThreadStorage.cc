// Mu2eG4PerThreadStorage holds the pointers to the UserActions and the various
// Manager classes, and other thread local objects for running G4 in MT mode.
//
// Author: Lisa Goodeough
// Date: 2019/07/09
//

#include "Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"

namespace mu2e {

  Mu2eG4PerThreadStorage::Mu2eG4PerThreadStorage(const Mu2eG4Config::Top& conf)
    : generatorModuleLabel{conf.generatorModuleLabel()}
    , multiStagePars{conf}
    , tvd{"", nullptr}
 {}

  //----------------------------------------------------------------
  void Mu2eG4PerThreadStorage::
  initializeEventInfo(art::Event* evt) {
    artEvent = evt;

    if(!(generatorModuleLabel == art::InputTag())) {
      artEvent->getByLabel(generatorModuleLabel, gensHandle);
    }

    // StepPointMCCollection of input hits from the previous simulation stage
    for(const auto& i : multiStagePars.genInputHits()) {
      genInputHits.emplace_back(evt->getValidHandle<StepPointMCCollection>(i));
    }

    // ProductID and ProductGetter for the SimParticleCollection.
    art::ProductID simPartId(evt->getProductID<SimParticleCollection>());
    art::EDProductGetter const* simProductGetter = evt->productGetter(simPartId);

    simParticleHelper = SimParticleHelper(multiStagePars.simParticleNumberOffset(), simPartId, evt, simProductGetter);
    simParticlePrimaryHelper = SimParticlePrimaryHelper(evt, simPartId, simProductGetter);

    if ( !gensHandle.isValid() && genInputHits.empty() ) {
      throw cet::exception("CONFIG")
        << "Error in PerThreadStorage::initializeEventInfo.  You are trying to run a G4job without an input for G4.\n";
    }

    if ( !gensHandle.isValid() && genInputHits.empty() ) {
      throw cet::exception("CONFIG")
        << "Error in PerThreadStorage::initializeEventInfo.  You are trying to run a G4job without an input for G4.\n";
    }
  }

  //----------------------------------------------------------------
  void Mu2eG4PerThreadStorage::
  putSensitiveDetectorData(art::EDProductGetter const* sim_product_getter) {

    std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> > steps_map =
      std::move(sensitiveDetectorSteps);

    for (std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> >::iterator i = steps_map.begin();
         i != steps_map.end(); ++i) {

      for ( StepPointMCCollection::iterator j=i->second->begin(); j!=i->second->end(); ++j ){

        StepPointMC& step = *j;
        step.simParticle() = art::Ptr<SimParticle>(step.simParticle().id(),
                                                   step.simParticle().key(),
                                                   sim_product_getter );
      }//for StepPointMCCollection::iterator

      artEvent->put(std::move(i->second), i->first);
    }//for (std::unordered_map...
  }

  //----------------------------------------------------------------
  void Mu2eG4PerThreadStorage::
  putCutsData(art::EDProductGetter const* sim_product_getter) {

    std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> > cuts_map = std::move(cutsSteps);

    for (std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> >::iterator i = cuts_map.begin();
         i != cuts_map.end(); ++i) {

      for ( StepPointMCCollection::iterator j=i->second->begin(); j!=i->second->end(); ++j ){
        StepPointMC& step = *j;

        step.simParticle() = art::Ptr<SimParticle>(step.simParticle().id(),
                                                   step.simParticle().key(),
                                                   sim_product_getter );
      }//for StepPointMCCollection::iterator

      artEvent->put(std::move(i->second), i->first);
    }
  }

  //----------------------------------------------------------------
  void Mu2eG4PerThreadStorage::clearData() {

    artEvent = nullptr;
    simParticleHelper = SimParticleHelper();
    simParticlePrimaryHelper = SimParticlePrimaryHelper();
    genInputHits.clear();
    gensHandle.clear();

    statG4 = nullptr;
    simPartCollection = nullptr;
    tvd.first = "";
    tvd.second = nullptr;
    mcTrajectories = nullptr;
    simRemapping = nullptr;
    extMonFNALHits = nullptr;
    sensitiveDetectorSteps.clear();
    cutsSteps.clear();
  }

  //----------------------------------------------------------------

} // end namespace mu2e
