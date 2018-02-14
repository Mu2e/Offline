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

#include "G4Threading.hh"


using namespace std;

namespace mu2e {

Mu2eG4RunAction::Mu2eG4RunAction(const bool using_MT,
                                 CLHEP::Hep3Vector const& origin_in_world,
                                 PhysicalVolumeHelper* phys_volume_helper,
                                 PhysicsProcessInfo* phys_process_info,
                                 TrackingAction *tracking_action,
                                 Mu2eG4SteppingAction *stepping_action,
                                 SensitiveDetectorHelper* sensitive_detectorhelper)
    :
    G4UserRunAction(),
    use_G4MT_(using_MT),
    originInWorld(origin_in_world),
    _physVolHelper(phys_volume_helper),
    _processInfo(phys_process_info),
    _trackingAction(tracking_action),
    _steppingAction(stepping_action),
    _sensitiveDetectorHelper(sensitive_detectorhelper)

    {
        //std::cout << "AT Mu2eG4RunAction c'tor" << std::endl;
    }
    
    
Mu2eG4RunAction::~Mu2eG4RunAction()
    {
        //std::cout << "AT Mu2eG4RunAction destructor" << std::endl;
    }
    
    
void Mu2eG4RunAction::BeginOfRunAction(const G4Run* aRun)
    {
        if (use_G4MT_ == true)//MT mode
        {
            
            //if (G4Threading::G4GetThreadId() == 0) {
            //    std::cout << "AT Mu2eG4RunAction::BeginOfRunAction in MT mode in Thread#0" << std::endl;
            //}
            

            if(IsMaster() == false)//do only in Worker Threads if in MT mode
            {
                _sensitiveDetectorHelper->registerSensitiveDetectors();
                /* THIS NEEDS TO BE ACTIVATED ONCE EVERYTHING IS WORKING
                 if (standardMu2eDetector_) _extMonFNALPixelSD =
                 dynamic_cast<ExtMonFNALPixelSD*>(G4SDManager::GetSDMpointer()
                 ->FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL()));
                 */


                _trackingAction->beginRun( _physVolHelper, _processInfo, originInWorld );
                _steppingAction->beginRun( _processInfo, originInWorld );
                
                //make this done only once per job/run according to ncalls in Mu2eG4_module?
                _steppingAction->finishConstruction();//once per thread
            }//if !Master
            
        }//MT mode
        else//sequential mode
        {
            //std::cout << "AT Mu2eG4RunAction::BeginOfRunAction in sequential mode" << std::endl;

            _sensitiveDetectorHelper->registerSensitiveDetectors();
            /* THIS NEEDS TO BE ACTIVATED ONCE EVERYTHING IS WORKING
             if (standardMu2eDetector_) _extMonFNALPixelSD =
             dynamic_cast<ExtMonFNALPixelSD*>(G4SDManager::GetSDMpointer()
             ->FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL()));
             */


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
        if (use_G4MT_ == false)//sequential mode
        {
        //std::cout << "calling Mu2eG4RA::EndOfRunAction()" << std::endl;
        _processInfo->endRun();
        // _physVolHelper.endRun(); is this really an endJob action that can be done in module?
        }
        
    }

    

}  // end namespace mu2e



