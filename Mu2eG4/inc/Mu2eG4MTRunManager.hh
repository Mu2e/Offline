#ifndef Mu2eG4_MTRunManager_hh
#define Mu2eG4_MTRunManager_hh
//
// Override the G4MTRunManager class so that the Mu2eG4 framework can drive
// the event loop.
//
// Original author Lisa Goodenough
//


// Included from Geant4
#include "Geant4/G4MTRunManager.hh"

//Mu2e includes
#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Offline/Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Offline/Mu2eG4/inc/SensitiveDetectorHelper.hh"

namespace G4 {class G4VUserPhysicsList;}

namespace mu2e {

  class Mu2eG4MasterRunAction;


  class Mu2eG4MTRunManager : public G4MTRunManager{

  public:

    Mu2eG4MTRunManager(const Mu2eG4Config::Top& conf, const Mu2eG4ResourceLimits& lim);
    virtual ~Mu2eG4MTRunManager();

    //we need our own versions of these functions in order to correctly control the event loop
    void initializeG4(int art_runnumber);
    void initializeKernelAndRM();
    void declarePhysicsAndGeometry();
    void initializeMasterRunAction();
    void stopG4();
    void terminateRun();
    G4bool SetUpEvent();

    inline G4VUserPhysicsList* getMasterPhysicsList() {return physicsList_;}

    inline void setPhysVolumeHelper(PhysicalVolumeHelper* phys_volume_helper) {physVolHelper_ = phys_volume_helper;}
    inline PhysicalVolumeHelper* getPhysVolumeHelper() {return physVolHelper_;}


  private:

    // Private and unimplemented to prevent copying.
    explicit Mu2eG4MTRunManager( Mu2eG4MTRunManager const & ) =  delete;
    explicit Mu2eG4MTRunManager( Mu2eG4MTRunManager && ) =  delete;
    Mu2eG4MTRunManager& operator=( Mu2eG4MTRunManager const & ) = delete;
    Mu2eG4MTRunManager& operator=( Mu2eG4MTRunManager && ) = delete;

    Mu2eG4Config::Top const conf_;
    Mu2eG4ResourceLimits const & mu2elimits_;

    bool m_managerInitialized;
    bool m_runTerminated;

    PhysicalVolumeHelper* physVolHelper_ = nullptr;
    SensitiveDetectorHelper sensitiveDetectorHelper_;
    Mu2eG4MasterRunAction* masterRunAction_ = nullptr;
    G4VUserPhysicsList* physicsList_ = nullptr;

    int rmvlevel_;
  };

} // end namespace mu2e

#endif /* Mu2eG4_MTRunManager_hh */
