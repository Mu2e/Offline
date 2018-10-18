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
#include "Mu2eG4/inc/ExtMonFNALPixelSD.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

//G4 includes
#include "G4Threading.hh"

//Mu2e includes
#include "GeometryService/inc/GeometryService.hh"


using namespace std;

namespace mu2e {

Mu2eG4RunAction::Mu2eG4RunAction(const bool using_MT,
                                 CLHEP::Hep3Vector const& origin_in_world,
                                 PhysicalVolumeHelper* phys_volume_helper,
                                 PhysicsProcessInfo* phys_process_info,
                                 TrackingAction *tracking_action,
                                 Mu2eG4SteppingAction *stepping_action,
                                 SensitiveDetectorHelper* sensitive_detectorhelper,
                                 ExtMonFNALPixelSD* extmon_FNAL_pixelSD)
    :
    G4UserRunAction(),
    use_G4MT_(using_MT),
    originInWorld(origin_in_world),
    _physVolHelper(phys_volume_helper),
    _processInfo(phys_process_info),
    _trackingAction(tracking_action),
    _steppingAction(stepping_action),
    _sensitiveDetectorHelper(sensitive_detectorhelper),
    _extMonFNALPixelSD(extmon_FNAL_pixelSD),
    standardMu2eDetector_((art::ServiceHandle<GeometryService>())->isStandardMu2eDetector())
    {}
    
    
Mu2eG4RunAction::~Mu2eG4RunAction()
    {}
    
    
void Mu2eG4RunAction::BeginOfRunAction(const G4Run* aRun)
    {
        if (use_G4MT_ == true){//MT mode
        
            if(IsMaster() == false) {//do only in Worker Threads if in MT mode
                
                _sensitiveDetectorHelper->registerSensitiveDetectors();
                
                _extMonFNALPixelSD = ( standardMu2eDetector_ &&
                                      _sensitiveDetectorHelper->extMonPixelsEnabled()) ?
                dynamic_cast<ExtMonFNALPixelSD*>(G4SDManager::GetSDMpointer()->
                                                 FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL()))
                : nullptr;
                
                _trackingAction->beginRun( _physVolHelper, _processInfo, originInWorld );
                _steppingAction->beginRun( _processInfo, originInWorld );
                
                //make this done only once per job/run according to ncalls in Mu2eG4_module?
                _steppingAction->finishConstruction();//once per thread
            }//if !Master
            
        }//MT mode
        else{//sequential mode

            _sensitiveDetectorHelper->registerSensitiveDetectors();
            
            _extMonFNALPixelSD = ( standardMu2eDetector_ &&
                                  _sensitiveDetectorHelper->extMonPixelsEnabled()) ?
            dynamic_cast<ExtMonFNALPixelSD*>(G4SDManager::GetSDMpointer()->
                                             FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL()))
            : nullptr;

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



