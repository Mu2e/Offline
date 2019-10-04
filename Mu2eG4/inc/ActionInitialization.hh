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
#include "Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"


//C++ includes
#include <vector>
#include <memory>

namespace fhicl { class ParameterSet; }

namespace CLHEP { class Hep3Vector; }

namespace mu2e {

    class SensitiveDetectorHelper;
    class IMu2eG4Cut;
    class PrimaryGeneratorAction;
    class GenEventBroker;
    class PerEventObjectsManager;
    class PhysicalVolumeHelper;


class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization(const fhicl::ParameterSet& pset,
                         std::vector< SensitiveDetectorHelper> &sensitive_detectorhelper_vector,
                         GenEventBroker* gen_eventbroker,
                         PhysicalVolumeHelper* phys_volume_helper,
                         const bool using_MT,
                         const int num_threads,
                         CLHEP::Hep3Vector const& origin_in_world,
                         Mu2eG4ResourceLimits const& mu2e_limits,
                         unsigned stage_offset_for_tracking_action
                         );

    virtual ~ActionInitialization();

    //BuildForMaster should be used for defining only the UserRunAction for the master thread.
    virtual void BuildForMaster() const;

    //Build should be used for defining user action classes for worker threads as well as for the sequential mode.
    virtual void Build() const;

    virtual G4VSteppingVerbose* InitializeSteppingVerbose() const;


   private:

    fhicl::ParameterSet const& pset_;

    //these are set using pset
    Mu2eG4TrajectoryControl trajectoryControl_;
    SimParticleCollectionPrinter simParticlePrinter_;
    std::vector<double> timeVDtimes_;
    Mu2eG4ResourceLimits mu2elimits_;


    std::vector < std::unique_ptr<IMu2eG4Cut> > stackingCutsVector;
    std::vector < std::unique_ptr<IMu2eG4Cut> > steppingCutsVector;
    std::vector < std::unique_ptr<IMu2eG4Cut> > commonCutsVector;

    GenEventBroker* _genEventBroker;
    PhysicalVolumeHelper* _physVolHelper;
    mutable PhysicsProcessInfo processInfo;

    const bool use_G4MT_;
    const int numthreads;
    CLHEP::Hep3Vector const& originInWorld;
    Mu2eG4ResourceLimits const& mu2eLimits;
    unsigned stageOffset;

    mutable std::vector< PerEventObjectsManager > perEvtObjManagerVector;
    mutable std::vector< PhysicsProcessInfo > physicsProcessInfoVector;
    std::vector< SensitiveDetectorHelper > &sensitiveDetectorHelperVector;

};

}  // end namespace mu2e
#endif /* Mu2eG4_ActionInitialization_hh */
