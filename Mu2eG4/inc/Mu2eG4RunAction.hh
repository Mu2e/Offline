#ifndef Mu2eG4_RunAction_hh
#define Mu2eG4_RunAction_hh
//
// Mu2eG4RunAction.hh provides declarations for the built-in RunAction
// class for the Mu2e G4 simulation.
//
// Author: Lisa Goodeough
// Date: 2017/08/07
//
//


//G4 includes
#include "Geant4/G4UserRunAction.hh"
#include "Geant4/G4Run.hh"
#include "Geant4/G4SDManager.hh"

#include "Mu2eG4/inc/Mu2eG4Config.hh"

namespace CLHEP { class Hep3Vector; }

namespace mu2e {

  class PhysicalVolumeHelper;
  class PhysicsProcessInfo;
  class Mu2eG4TrackingAction;
  class Mu2eG4SteppingAction;
  class SensitiveDetectorHelper;

  class Mu2eG4RunAction : public G4UserRunAction
  {
  public:
    Mu2eG4RunAction(const Mu2eG4Config::Debug& debug,
                    CLHEP::Hep3Vector const&,
                    PhysicalVolumeHelper*,
                    PhysicsProcessInfo*,
                    Mu2eG4TrackingAction*,
                    Mu2eG4SteppingAction*,
                    SensitiveDetectorHelper*
                    );

    virtual ~Mu2eG4RunAction();

    //methods
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

  private:
    //data members

    Mu2eG4Config::Debug debug_;
    CLHEP::Hep3Vector const& originInWorld;

    PhysicalVolumeHelper* _physVolHelper;
    PhysicsProcessInfo* _processInfo;
    Mu2eG4TrackingAction* _trackingAction;
    Mu2eG4SteppingAction* _steppingAction;

    SensitiveDetectorHelper* _sensitiveDetectorHelper;

  };

}  // end namespace mu2e
#endif /* Mu2eG4_RunAction_hh */
