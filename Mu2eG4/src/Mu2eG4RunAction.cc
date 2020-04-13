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

//G4 includes
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;

namespace mu2e {

  Mu2eG4RunAction::Mu2eG4RunAction(const Mu2eG4Config::Debug& debug,
                                   CLHEP::Hep3Vector const& origin_in_world,
                                   PhysicalVolumeHelper* phys_volume_helper,
                                   PhysicsProcessInfo* phys_process_info,
                                   TrackingAction *tracking_action,
                                   Mu2eG4SteppingAction *stepping_action,
                                   SensitiveDetectorHelper* sensitive_detectorhelper
                                   )
  :
  G4UserRunAction(),
    debug_(debug),
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

    if (debug_.diagLevel() > 0) {
      G4cout << "Mu2eG4RunAction " << __func__ << " : G4Run: " << aRun->GetRunID() << G4endl;
    }

    // run managers are thread local
    G4RunManagerKernel const * rmk = G4RunManagerKernel::GetRunManagerKernel();
    G4TrackingManager* tm  = rmk->GetTrackingManager();
    tm->SetVerboseLevel(debug_.trackingVerbosityLevel());
    G4SteppingManager* sm  = tm->GetSteppingManager();
    sm->SetVerboseLevel(debug_.steppingVerbosityLevel());

    _sensitiveDetectorHelper->registerSensitiveDetectors();

    if (!_physVolHelper->helperIsInitialized())
      {
        _physVolHelper->beginRun();//map w/~20,000 entries
      }

    _processInfo->beginRun();

    _trackingAction->beginRun( _physVolHelper, _processInfo, originInWorld );
    _steppingAction->beginRun( _processInfo, originInWorld );
    _steppingAction->finishConstruction();
    G4Navigator* navigator =
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    navigator->CheckMode(debug_.navigatorCheckMode());
    navigator->SetVerboseLevel(debug_.navigatorVerbosityLevel());

  }


  void Mu2eG4RunAction::EndOfRunAction(const G4Run* aRun)
  {
    _processInfo->endRun();
  }

}  // end namespace mu2e
