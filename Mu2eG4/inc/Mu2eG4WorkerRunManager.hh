#ifndef Mu2eG4_WorkerRunManager_hh
#define Mu2eG4_WorkerRunManager_hh
//
// Override the G4WorkerRunManager class so that the Mu2eG4 framework can drive
// the event loop.
//
// Original author Lisa Goodenough
//


// Included from Geant4
#include "Geant4/G4WorkerRunManager.hh"

//Mu2e includes
#include "Mu2eG4/inc/Mu2eG4Config.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"

// C++ includes
#include <thread>


namespace art { class Event; }
namespace CLHEP { class Hep3Vector; }

namespace mu2e {

  class Mu2eG4PerThreadStorage;

  class Mu2eG4MTRunManager;
  class PrimaryGeneratorAction;
  class Mu2eG4SteppingAction;
  class Mu2eG4TrackingAction;
  class ExtMonFNALPixelSD;
  class Mu2eG4IOConfigHelper;

  class Mu2eG4WorkerRunManager : public G4WorkerRunManager{

  public:

    Mu2eG4WorkerRunManager(const Mu2eG4Config::Top& conf, const Mu2eG4IOConfigHelper& ioconf, std::thread::id worker_ID);
    virtual ~Mu2eG4WorkerRunManager();

    //**********************************************************
    //These functions help us control the event loop
    void initializeThread(Mu2eG4MTRunManager* mRM,
                          const G4ThreeVector& origin_in_world);
    void initializeUserActions(const G4ThreeVector& origin_in_world);
    void initializeRun(art::Event* const art_event);
    void processEvent(const art::EventID&);
    G4Event* generateEvt(const art::EventID&);

    inline bool workerRMInitialized() const { return m_managerInitialized; }

    Mu2eG4PerThreadStorage* getMu2eG4PerThreadStorage() {
      return perThreadObjects_.get();
    }

  private:

    Mu2eG4Config::Top conf_;

    bool m_managerInitialized;
    bool m_steppingVerbose;
    int m_mtDebugOutput;
    int rmvlevel_;
    std::string salt_;

    std::unique_ptr<Mu2eG4PerThreadStorage> perThreadObjects_;

    Mu2eG4MTRunManager* masterRM;

    std::thread::id workerID_;

    PhysicsProcessInfo physicsProcessInfo_;
    SensitiveDetectorHelper sensitiveDetectorHelper_;
    ExtMonFNALPixelSD* extMonFNALPixelSD_;

    PrimaryGeneratorAction* genAction_;
    Mu2eG4SteppingAction* steppingAction_;
    Mu2eG4TrackingAction* trackingAction_;

  };

} // end namespace mu2e

#endif /* Mu2eG4_WorkerRunManager_hh */
