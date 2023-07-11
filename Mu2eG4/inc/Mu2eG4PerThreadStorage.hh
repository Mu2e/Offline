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
#include "Offline/Mu2eG4/inc/Mu2eG4IOConfigHelper.hh"
#include "Offline/Mu2eG4/inc/SimParticleHelper.hh"
#include "Offline/Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "Offline/Mu2eG4/inc/IMu2eG4Cut.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/StatusG4.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/SimParticleRemapping.hh"
#include "Offline/MCDataProducts/inc/ExtMonFNALSimHit.hh"


// C++ includes
#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include <optional>

//art includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Utilities/InputTag.h"

namespace art { class EDProductGetter; }

namespace mu2e {

  struct Mu2eG4PerThreadStorage
  {
    explicit Mu2eG4PerThreadStorage(const Mu2eG4IOConfigHelper& ioconf);

    void initializeEventInfo(art::Event* evt, unsigned simStage);

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

    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
    // functions to move the data into the art::Event

    bool eventPassed() const { return bool(statG4); }
    void putDataIntoEvent();

    void putStepPointMCCollections(std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> > &&steps_map);

    void putSensitiveDetectorData();

    void clearData();

    /////////////////////////////////////////////////////////////
    // Job-level configuration
    Mu2eG4IOConfigHelper ioconf;

    // run-level data members
    art::RunNumber_t currentRunNumber = 0;

    bool runTerminated = false;

    // event-level data members
    art::Event* artEvent = nullptr;
    // delay initialization of the helpers using std::optional
    std::optional<SimParticleHelper> simParticleHelper;
    std::optional<SimParticlePrimaryHelper> simParticlePrimaryHelper;

    // output data products
    std::unique_ptr<StatusG4> statG4{nullptr};
    std::unique_ptr<SimParticleCollection> simPartCollection = nullptr;
    std::unique_ptr<StepPointMCCollection> tvd_collection = nullptr;
    std::unique_ptr<MCTrajectoryCollection> mcTrajectories = nullptr;
    std::unique_ptr<SimParticleRemapping> simRemapping = nullptr;
    std::unique_ptr<ExtMonFNALSimHitCollection> extMonFNALHits = nullptr;

    std::unordered_map< std::string, std::unique_ptr<StepPointMCCollection> > sensitiveDetectorSteps;

    std::unique_ptr<IMu2eG4Cut> stackingCuts = nullptr;
    std::unique_ptr<IMu2eG4Cut> steppingCuts = nullptr;
    std::unique_ptr<IMu2eG4Cut> commonCuts = nullptr;
  };

}  // end namespace mu2e
#endif /* Mu2eG4_PerThreadStorage_hh */
