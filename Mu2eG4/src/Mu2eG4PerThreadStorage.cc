// Mu2eG4PerThreadStorage holds the pointers to the UserActions and the various
// Manager classes, and other thread local objects for running G4 in MT mode.
//
// Author: Lisa Goodeough
// Date: 2019/07/09
//

#include "Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"

namespace mu2e {

  Mu2eG4PerThreadStorage::Mu2eG4PerThreadStorage(const Mu2eG4Config::Top& conf)
    : generatorModuleLabel{conf.generatorModuleLabel()}
    , multiStagePars{conf}
    , timeVD_enabled(conf.SDConfig().TimeVD().enabled())
    , produceMCTrajectories{conf.TrajectoryControl().produce()}
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

    // Output collections
    simPartCollection = std::unique_ptr<SimParticleCollection>( new SimParticleCollection );

    if(timeVD_enabled) {
      tvd_collection = std::unique_ptr<StepPointMCCollection>( new StepPointMCCollection );
    }

    if(produceMCTrajectories) {
      mcTrajectories = std::unique_ptr<MCTrajectoryCollection>( new MCTrajectoryCollection );
    }

    if(multiStagePars.multiStage()) {
      simRemapping = std::unique_ptr<SimParticleRemapping>( new SimParticleRemapping );
    }

  }

  //----------------------------------------------------------------
  void Mu2eG4PerThreadStorage::
  putStepPointMCCollections(std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> > &&steps_map) {

    art::ProductID simPartId(artEvent->getProductID<SimParticleCollection>());
    art::EDProductGetter const* simProductGetter = artEvent->productGetter(simPartId);

    for (auto i = steps_map.begin(); i != steps_map.end(); ++i) {
      for (auto j=i->second->begin(); j!=i->second->end(); ++j ){

        StepPointMC& step = *j;
        step.simParticle() = art::Ptr<SimParticle>(step.simParticle().id(),
                                                   step.simParticle().key(),
                                                   simProductGetter);
      }

      artEvent->put(std::move(i->second), i->first);
    }
  }

  //----------------------------------------------------------------
  void Mu2eG4PerThreadStorage::putSensitiveDetectorData() {
    putStepPointMCCollections(std::move(sensitiveDetectorSteps));
  }

  //----------------------------------------------------------------
  void Mu2eG4PerThreadStorage::putCutsData() {
    putStepPointMCCollections(std::move(cutsSteps));
  }

  //----------------------------------------------------------------
  void Mu2eG4PerThreadStorage::putDataIntoEvent() {
      artEvent->put(std::move(statG4));
      artEvent->put(std::move(simPartCollection));

      putSensitiveDetectorData();
      putCutsData();

      if(timeVD_enabled) {
        static const StepInstanceName timeVD(StepInstanceName::timeVD);
        artEvent->put(std::move(tvd_collection), timeVD.name());
      }

      if(produceMCTrajectories) {
        artEvent->put(std::move(mcTrajectories));
      }

      if(multiStagePars.multiStage()) {
        artEvent->put(std::move(simRemapping));
      }

      if(extMonFNALHits) {
        artEvent->put(std::move(extMonFNALHits));
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
    tvd_collection = nullptr;
    mcTrajectories = nullptr;
    simRemapping = nullptr;
    extMonFNALHits = nullptr;
    sensitiveDetectorSteps.clear();
    cutsSteps.clear();
  }

  //----------------------------------------------------------------

} // end namespace mu2e
