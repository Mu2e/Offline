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
#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "G4SDManager.hh"

#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {
    
    
    class PhysicalVolumeHelper;
    class PhysicsProcessInfo;
    class TrackingAction;
    class Mu2eG4SteppingAction;
    class SensitiveDetectorHelper;
    class ExtMonFNALPixelSD;

class Mu2eG4RunAction : public G4UserRunAction
{
  public:
    Mu2eG4RunAction(const bool,
                    CLHEP::Hep3Vector const&,
                    PhysicalVolumeHelper*,
                    PhysicsProcessInfo*,
                    TrackingAction*,
                    Mu2eG4SteppingAction*,
                    SensitiveDetectorHelper*,
                    ExtMonFNALPixelSD*
                    );
    
    virtual ~Mu2eG4RunAction();
    
    //methods
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);
    
    
   private:
    //data members
    const bool use_G4MT_;
    CLHEP::Hep3Vector const& originInWorld;
    
    PhysicalVolumeHelper* _physVolHelper;
    PhysicsProcessInfo* _processInfo;
    TrackingAction* _trackingAction;
    Mu2eG4SteppingAction* _steppingAction;
    
    SensitiveDetectorHelper* _sensitiveDetectorHelper;
    ExtMonFNALPixelSD* _extMonFNALPixelSD;
    
    const bool standardMu2eDetector_;
    
    
};

}  // end namespace mu2e
#endif /* Mu2eG4_RunAction_hh */


