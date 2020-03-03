#ifndef Mu2eG4_WorkerRunManager_hh
#define Mu2eG4_WorkerRunManager_hh
//
// Override the G4WorkerRunManager class so that the Mu2eG4 framework can drive
// the event loop.
//
// Original author Lisa Goodenough
//


// Included from Geant4
#include "G4WorkerRunManager.hh"

//Mu2e includes
#include "Mu2eG4/inc/Mu2eG4Config.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eG4/inc/Mu2eG4MultiStageParameters.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"

// C++ includes
#include <thread>


namespace art { class Event; }
namespace CLHEP { class Hep3Vector; }

namespace mu2e {

  class Mu2eG4PerThreadStorage;

  class Mu2eG4MTRunManager;
  class IMu2eG4Cut;
  class PrimaryGeneratorAction;
  class Mu2eG4SteppingAction;
  class TrackingAction;
  class ExtMonFNALPixelSD;

  class Mu2eG4WorkerRunManager : public G4WorkerRunManager{

  public:

    Mu2eG4WorkerRunManager(const Mu2eG4Config::Top& conf, std::thread::id worker_ID);
    virtual ~Mu2eG4WorkerRunManager();

    //**********************************************************
    //These functions help us control the event loop
    void initializeThread(Mu2eG4MTRunManager* mRM,
                          const G4ThreeVector& origin_in_world);
    void initializeUserActions(const G4ThreeVector& origin_in_world);
    void initializeRun(art::Event* art_event);
    void processEvent(art::Event*);

    inline bool workerRMInitialized() const { return m_managerInitialized; }

    Mu2eG4PerThreadStorage* getMu2eG4PerThreadStorage() {
      return perThreadObjects_.get();
    }

  private:

    Mu2eG4Config::Top conf_;

    bool m_managerInitialized;
    bool m_steppingVerbose;
    bool m_mtDebugOutput;
    int rmvlevel_;

    std::unique_ptr<Mu2eG4PerThreadStorage> perThreadObjects_;

    Mu2eG4MTRunManager* masterRM;

    std::thread::id workerID_;
    Mu2eG4ResourceLimits mu2elimits_;
    Mu2eG4TrajectoryControl trajectoryControl_;
    Mu2eG4MultiStageParameters multiStagePars_;

    PhysicsProcessInfo physicsProcessInfo_;
    SensitiveDetectorHelper sensitiveDetectorHelper_;
    ExtMonFNALPixelSD* extMonFNALPixelSD_;
    std::unique_ptr<IMu2eG4Cut> stackingCuts_;
    std::unique_ptr<IMu2eG4Cut> steppingCuts_;
    std::unique_ptr<IMu2eG4Cut> commonCuts_;

    PrimaryGeneratorAction* genAction_;
    Mu2eG4SteppingAction* steppingAction_;
    TrackingAction* trackingAction_;

  };

} // end namespace mu2e

#endif /* Mu2eG4_WorkerRunManager_hh */
