//
//
// Author: Lisa Goodeough
// Date: 2017/08/07
//
//

//Mu2e includes
#include "Mu2eG4/inc/Mu2eG4RunAction.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/Mu2eG4SteppingAction.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

//art includes
#include "fhiclcpp/ParameterSet.h"

//G4 includes
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;

namespace mu2e {

Mu2eG4RunAction::Mu2eG4RunAction(const fhicl::ParameterSet& pset,
                                 const bool using_MT,
                                 CLHEP::Hep3Vector const& origin_in_world,
                                 PhysicalVolumeHelper* phys_volume_helper,
                                 PhysicsProcessInfo* phys_process_info,
                                 TrackingAction *tracking_action,
                                 Mu2eG4SteppingAction *stepping_action,
                                 SensitiveDetectorHelper* sensitive_detectorhelper
                                 )
    :
    G4UserRunAction(),
    pset_(pset),
    use_G4MT_(using_MT),
    originInWorld(origin_in_world),
    _physVolHelper(phys_volume_helper),
    _processInfo(phys_process_info),
    _trackingAction(tracking_action),
    _steppingAction(stepping_action),
    _sensitiveDetectorHelper(sensitive_detectorhelper)
    {}

Mu2eG4RunAction::~Mu2eG4RunAction()
    {}
    
    
void Mu2eG4RunAction::BeginOfRunAction(const G4Run* aRun)
    {

      if (pset_.get<int>("debug.diagLevel",0)>1) {
        G4cout << "Mu2eG4RunAction " << __func__ << " : G4Run: " << aRun->GetRunID() << G4endl;
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

        if (use_G4MT_ == true){//MT mode

          //BeginOfRunAction is called from Worker Threads only

          _sensitiveDetectorHelper->registerSensitiveDetectors();
            
          _trackingAction->beginRun( _physVolHelper, _processInfo, originInWorld );
          _steppingAction->beginRun( _processInfo, originInWorld );
                
          //make this done only once per job/run according to ncalls in Mu2eG4_module?
          _steppingAction->finishConstruction();//once per thread
            
        }//MT mode
        else{//sequential mode

            _sensitiveDetectorHelper->registerSensitiveDetectors();
            
            _physVolHelper->beginRun();//map w/~20,000 entries
            _processInfo->beginRun();
            
            _trackingAction->beginRun( _physVolHelper, _processInfo, originInWorld );
            _steppingAction->beginRun( _processInfo, originInWorld );
            
            //make this done only once per job/run according to ncalls in Mu2eG4_module?
            _steppingAction->finishConstruction();//once per thread
            
        }//sequential mode
        
    }
    
    
void Mu2eG4RunAction::EndOfRunAction(const G4Run* aRun)
    {
        if (use_G4MT_ == false){//sequential mode
            _processInfo->endRun();
            //_physVolHelper.endRun(); is this really an endJob action that can be done in module?
        }
    }

}  // end namespace mu2e



