#ifndef Mu2eG4_ActionInitialization_hh
#define Mu2eG4_ActionInitialization_hh
//
// ActionInitialization.hh provides declarations for the built-in action
// initialization for the Mu2e G4 simulation. In both BuildForMaster and Build,
// it instantiates user action classes and registers user actions through the
// SetUserAction method defined in the G4VUserActionInitialization base class.
//
// Author: Lisa Goodeough
// Date: 2017/05/18
//
//

//G4 includes
#include "G4VUserActionInitialization.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

//Mu2e includes
#include "Mu2eG4/inc/Mu2eG4Config.hh"
#include "Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"


//C++ includes
#include <vector>
#include <memory>

namespace CLHEP { class Hep3Vector; }

namespace mu2e {

  class SensitiveDetectorHelper;
  class IMu2eG4Cut;
  class PrimaryGeneratorAction;
  class Mu2eG4PerThreadStorage;
  class PhysicalVolumeHelper;


  class ActionInitialization : public G4VUserActionInitialization
  {
  public:
    ActionInitialization(const Mu2eG4Config::Top& conf,
                         SensitiveDetectorHelper* sensitive_detectorhelper,
                         Mu2eG4PerThreadStorage* per_thread_storage,
                         PhysicalVolumeHelper* phys_volume_helper,
                         CLHEP::Hep3Vector const& origin_in_world,
                         unsigned stage_offset_for_tracking_action
                         );

    virtual ~ActionInitialization();

    //BuildForMaster should be used for defining only the UserRunAction for the master thread.
    virtual void BuildForMaster() const;

    //Build should be used for defining user action classes for worker threads as well as for the sequential mode.
    virtual void Build() const;

    virtual G4VSteppingVerbose* InitializeSteppingVerbose() const;


  private:

    Mu2eG4Config::Top conf_;

    Mu2eG4TrajectoryControl trajectoryControl_;
    std::vector<double> timeVDtimes_;
    Mu2eG4ResourceLimits mu2eLimits_;

    std::unique_ptr<IMu2eG4Cut> stackingCuts_;
    std::unique_ptr<IMu2eG4Cut> steppingCuts_;
    std::unique_ptr<IMu2eG4Cut> commonCuts_;

    SensitiveDetectorHelper* sensitiveDetectorHelper_;
    Mu2eG4PerThreadStorage*  perThreadStorage_;
    PhysicalVolumeHelper* physVolHelper_;

    mutable PhysicsProcessInfo physicsProcessInfo_;

    CLHEP::Hep3Vector const& originInWorld_;
    unsigned stageOffset_;
  };

}  // end namespace mu2e
#endif /* Mu2eG4_ActionInitialization_hh */
