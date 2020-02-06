//
//
// Author: Lisa Goodeough
// Date: 2017/08/07
//
//

//G4 includes
#include "G4VUserPhysicsList.hh"
#include "G4RunManager.hh"

//Mu2e includes
#include "Mu2eG4/inc/Mu2eG4MasterRunAction.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"

using namespace std;

namespace mu2e {

  Mu2eG4MasterRunAction::Mu2eG4MasterRunAction(int diagLevel,
                                               PhysicalVolumeHelper* phys_volume_helper)
    :
    G4UserRunAction(),
    diagLevel_(diagLevel),
    physVolHelper_(phys_volume_helper)
  {}


  Mu2eG4MasterRunAction::~Mu2eG4MasterRunAction()
  {}


  void Mu2eG4MasterRunAction::BeginOfRunAction(const G4Run* aRun) {
    //this can be called in ActionInitialization::BuildForMaster()
    //however, for our version of MT running, we are not utilizing the ActionInitialization class
  }


  void Mu2eG4MasterRunAction::MasterBeginRunAction() {

    if (diagLevel_ > 0) {
      G4cout << "Mu2eG4MasterRunAction " << __func__ << " called " << G4endl;
    }

    //this class is ONLY called in MT mode
    //we want these actions performed only in the Master thread
    physVolHelper_->beginRun();//map w/~20,000 entries

  }


  void Mu2eG4MasterRunAction::MasterEndRunAction()
  {}


}  // end namespace mu2e
