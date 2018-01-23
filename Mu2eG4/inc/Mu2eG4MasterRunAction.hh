#ifndef Mu2eG4_MasterRunAction_hh
#define Mu2eG4_MasterRunAction_hh
//
// Mu2eG4MasterRunAction.hh provides declarations for the built-in RunAction
// class for the Mu2e G4 simulation for the Master Thread in MT mode.
// This class is instantiated only in ActionInitialization::BuildForMaster()
//
// Author: Lisa Goodeough
// Date: 2017/08/07
//
//


//G4 includes
#include "G4UserRunAction.hh"
#include "G4Run.hh"

#include <vector>


namespace mu2e {
    
    class PhysicalVolumeHelper;
    class PhysicsProcessInfo;

class Mu2eG4MasterRunAction : public G4UserRunAction
{
  public:
    //Mu2eG4MasterRunAction(PhysicalVolumeHelper*, PhysicsProcessInfo*);

    Mu2eG4MasterRunAction(PhysicalVolumeHelper*,
                          std::vector< PhysicsProcessInfo >*);

    virtual ~Mu2eG4MasterRunAction();
    
    //methods
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);
    
    
   private:
    //data members
    
    PhysicalVolumeHelper* _physVolHelper;
    //PhysicsProcessInfo* _processInfo;
    std::vector< PhysicsProcessInfo >* PhysicsProcessInfoVector;
    
};

}  // end namespace mu2e
#endif /* Mu2eG4_MasterRunAction_hh */


