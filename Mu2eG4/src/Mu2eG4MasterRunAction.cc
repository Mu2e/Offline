//
//
// Author: Lisa Goodeough
// Date: 2017/08/07
//
//

//G4 includes
#include "G4VUserPhysicsList.hh"
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"

//art includes
#include "fhiclcpp/ParameterSet.h"

//Mu2e includes
#include "Mu2eG4/inc/Mu2eG4MasterRunAction.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"

using namespace std;

namespace mu2e {

Mu2eG4MasterRunAction::Mu2eG4MasterRunAction(const fhicl::ParameterSet& pset,
                                             PhysicalVolumeHelper* phys_volume_helper,
                                             std::vector< PhysicsProcessInfo >* phys_process_info_vec
                                             )
    :
    G4UserRunAction(),
    pset_(pset),
    _physVolHelper(phys_volume_helper),
    PhysicsProcessInfoVector(phys_process_info_vec)
    {}


Mu2eG4MasterRunAction::~Mu2eG4MasterRunAction()
    {}


void Mu2eG4MasterRunAction::BeginOfRunAction(const G4Run* aRun)
    {

      if (pset_.get<int>("debug.diagLevel",0)>0) {
        G4cout << "Mu2eG4MasterRunAction " << __func__ << " called " << G4endl;
      }

      // run managers are thread local
      G4RunManagerKernel const * rmk = G4RunManagerKernel::GetRunManagerKernel();
      G4TrackingManager* tm  = rmk->GetTrackingManager();
      tm->SetVerboseLevel(pset_.get<int>("debug.trackingVerbosityLevel",0));
      G4SteppingManager* sm  = tm->GetSteppingManager();
      sm->SetVerboseLevel(pset_.get<int>("debug.steppingVerbosityLevel",0));
      G4Navigator* navigator =
	G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
      navigator->CheckMode(pset_.get<bool>("debug.navigatorCheckMode",false));

        //this class is ONLY called in MT mode
        //we want these actions performed only in the Master thread
        _physVolHelper->beginRun();//map w/~20,000 entries
        
        for (unsigned i = 0; i < PhysicsProcessInfoVector->size(); i++) {
            PhysicsProcessInfoVector->at(i).beginRun();
        }

    }
    
    
void Mu2eG4MasterRunAction::EndOfRunAction(const G4Run* aRun)
    {        
        for (unsigned i = 0; i < PhysicsProcessInfoVector->size(); i++) {
            PhysicsProcessInfoVector->at(i).endRun();
        }
    }

}  // end namespace mu2e



