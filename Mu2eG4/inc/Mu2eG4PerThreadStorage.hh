#ifndef Mu2eG4_PerThreadStorage_hh
#define Mu2eG4_PerThreadStorage_hh
//
// Mu2eG4PerThreadStorage.hh holds the pointers to the UserActions and the various
// Manager classes, and other thread local objects for running G4 in MT mode.
//
// Author: Lisa Goodeough
// Date: 2019/07/09
//

//Mu2e includes
#include "Mu2eG4/inc/Mu2eG4Config.hh"
#include "Mu2eG4/inc/Mu2eG4MultiStageParameters.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"

// C++ includes
#include <memory>
#include <vector>
#include <string>
#include <unordered_map>

//art includes
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Utilities/InputTag.h"

namespace art { class EDProductGetter; }

namespace mu2e {

  typedef std::vector< art::ValidHandle<StepPointMCCollection> > HitHandles;

  struct Mu2eG4PerThreadStorage
  {
    explicit Mu2eG4PerThreadStorage(const Mu2eG4Config::Top& conf);

    void initializeEventInfo(art::Event* evt);

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    // functions to get the event data from the EventAction
    void insertMCTrajectoryCollection(std::unique_ptr<MCTrajectoryCollection> mc_trajectories) {
      mcTrajectories = std::move(mc_trajectories);
    }

    void insertSDStepPointMC(std::unique_ptr<StepPointMCCollection> step_point_mc,
                             std::string instance_name) {
      sensitiveDetectorSteps[instance_name] = std::move(step_point_mc);

    }

    void insertCutsStepPointMC(std::unique_ptr<StepPointMCCollection> step_point_mc,
                               std::string instance_name) {
      cutsSteps[instance_name] = std::move(step_point_mc);
    }

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    // functions to move the data into the art::Event

    bool eventPassed() const { return bool(statG4); }
    void putDataIntoEvent();

    void putStepPointMCCollections(std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> > &&steps_map);

    void putSensitiveDetectorData();

    void putCutsData();

    void clearData();

    /////////////////////////////////////////////////////////////
    const art::InputTag generatorModuleLabel;
    const Mu2eG4MultiStageParameters multiStagePars;
    const bool timeVD_enabled;
    const bool produceMCTrajectories;

    // run-level data members
    art::RunNumber_t currentRunNumber = 0;

    bool runTerminated = false;

    // event-level data members
    art::Event* artEvent = nullptr;
    SimParticleHelper simParticleHelper;
    SimParticlePrimaryHelper simParticlePrimaryHelper;
    HitHandles genInputHits;
    art::Handle<GenParticleCollection> gensHandle;


    // output data products
    std::unique_ptr<StatusG4> statG4{nullptr};
    std::unique_ptr<SimParticleCollection> simPartCollection = nullptr;
    std::unique_ptr<StepPointMCCollection> tvd_collection;
    std::unique_ptr<MCTrajectoryCollection> mcTrajectories = nullptr;
    std::unique_ptr<SimParticleRemapping> simRemapping = nullptr;
    std::unique_ptr<ExtMonFNALSimHitCollection> extMonFNALHits = nullptr;

    std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> > sensitiveDetectorSteps;
    std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> > cutsSteps;
  };

}  // end namespace mu2e
#endif /* Mu2eG4_PerThreadStorage_hh */
