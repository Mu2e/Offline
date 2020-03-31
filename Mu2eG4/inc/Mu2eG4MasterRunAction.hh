#ifndef Mu2eG4_MasterRunAction_hh
#define Mu2eG4_MasterRunAction_hh
//
// Mu2eG4MasterRunAction.hh provides declarations for the built-in RunAction
// class for the Mu2e G4 simulation for the Master Thread in MT mode.
// This class is instantiated only in ActionInitialization::BuildForMaster()
//
// Author: Lisa Goodenough
// Date: 2017/08/07
//
//


//G4 includes
#include "G4UserRunAction.hh"
#include "G4Run.hh"

#include <vector>

namespace mu2e {

  class PhysicalVolumeHelper;

  class Mu2eG4MasterRunAction : public G4UserRunAction
  {

  public:
    Mu2eG4MasterRunAction(int diagLevel_,
                          PhysicalVolumeHelper*);

    virtual ~Mu2eG4MasterRunAction();

    //methods
    virtual void BeginOfRunAction(const G4Run* aRun);

    void MasterBeginRunAction();
    void MasterEndRunAction();

  private:
    //data members
    int diagLevel_;
    PhysicalVolumeHelper* physVolHelper_;

  };

}  // end namespace mu2e
#endif /* Mu2eG4_MasterRunAction_hh */
