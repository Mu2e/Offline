//
// ActionInitialization.cc provides implementation of Mu2e G4's built-in action initialization.
//
// Author: Lisa Goodenough
// Date: 2017/05/08
//
//

//Mu2e includes
#include "Mu2eG4/inc/ActionInitialization.hh"
#include "Mu2eG4/inc/PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/Mu2eG4StackingAction.hh"
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/Mu2eG4SteppingAction.hh"
#include "Mu2eG4/inc/Mu2eG4EventAction.hh"
#include "Mu2eG4/inc/Mu2eG4RunAction.hh"
#include "Mu2eG4/inc/Mu2eG4MasterRunAction.hh"
#include "Mu2eG4/inc/ExtMonFNALPixelSD.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/IMu2eG4Cut.hh"
#include "Mu2eG4/inc/GenEventBroker.hh"
#include "Mu2eG4/inc/PerEventObjectsManager.hh"
#include "Mu2eG4/inc/SteppingVerbose.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"

//C++ includes
#include <iostream>

//G4 includes
//#include "G4Threading.hh"


namespace mu2e {
    
    ActionInitialization::ActionInitialization(const fhicl::ParameterSet& pset,
                                               ExtMonFNALPixelSD* extmon_FNAL_pixelSD,
                                               std::vector< SensitiveDetectorHelper> &sensitive_detectorhelper_vector,
                                               IMu2eG4Cut& stacking_cuts,
                                               IMu2eG4Cut& stepping_cuts,
                                               IMu2eG4Cut& common_cuts,
                                               GenEventBroker* gen_eventbroker,
                                               PhysicalVolumeHelper* phys_volume_helper,
                                               const bool using_MT,
                                               const int num_threads,
                                               CLHEP::Hep3Vector const& origin_in_world
                                               )
        :
        G4VUserActionInitialization(),
        pset_(pset),
        trajectoryControl_(pset.get<fhicl::ParameterSet>("TrajectoryControl")),
        simParticlePrinter_(pset.get<fhicl::ParameterSet>("SimParticlePrinter", SimParticleCollectionPrinter::defaultPSet())),
        timeVDtimes_(pset.get<std::vector<double> >("SDConfig.TimeVD.times")),
        mu2elimits_(pset.get<fhicl::ParameterSet>("ResourceLimits")),
        _extMonFNALPixelSD(extmon_FNAL_pixelSD),
        _stackingCuts(stacking_cuts),
        _steppingCuts(stepping_cuts),
        _commonCuts(common_cuts),
        _genEventBroker(gen_eventbroker),
        _physVolHelper(phys_volume_helper),
        processInfo(),
        use_G4MT_(using_MT),
        numthreads(num_threads),
        originInWorld(origin_in_world),
        sensitiveDetectorHelperVector(sensitive_detectorhelper_vector)
        {
            
            //_sensitiveDetectorHelper = &(sensitiveDetectorHelperVector.at(0));
            
            perEvtObjManagerVector.reserve(numthreads);
            physicsProcessInfoVector.reserve(numthreads);
            
            //std::cout << "size of perEvtObjManagerVector: " << perEvtObjManagerVector.size() << std::endl;
            for (int i = 0; i < numthreads ; i++) {
                //std::cout << "calling perEvtObjManagerVector.emplace_back for the " << i <<  "th time" << std::endl;
                //perEvtObjManagerVector.emplace_back(PerEventObjectsManager(pset_));
                perEvtObjManagerVector.emplace_back(pset_,numthreads);
            }
            //std::cout << "size of perEvtObjManagerVector after emplace_back: " << perEvtObjManagerVector.size() << std::endl;

            //std::cout << "Address of this PEOM vector is " << &perEvtObjManagerVector << std::endl;
            //change this to # G4events? AI constructor called only once in both MT & seq mode
            
            
            for (int i = 0; i < numthreads ; i++) {
                physicsProcessInfoVector.emplace_back();
            }

            
        }
    
    
    ActionInitialization::~ActionInitialization()
    {        
        std::cout << "AT ActionInitialization destructor" << std::endl;
        
    }

    
// used for defining only the UserRunAction for the master thread
    void ActionInitialization::BuildForMaster() const
    {
        SetUserAction( new Mu2eG4MasterRunAction(_physVolHelper, &physicsProcessInfoVector) );
    }
    
    
// used for defining user action classes for worker threads as well as for the sequential mode.
    void ActionInitialization::Build() const
    {
//        std::cout << "WE ARE INSIDE BUILD()! in the Master Thread: " << G4Threading::IsMasterThread()
//        << ", from THREAD #" << G4Threading::G4GetThreadId() << std::endl;
        
        if (G4Threading::IsMasterThread() == true){
            std::cout << "WE ARE INSIDE BUILD() in the Master Thread" << std::endl;
        }
        else{
            //std::cout << "WE ARE INSIDE BUILD() from Thread #" << G4Threading::G4GetThreadId() << std::endl;
        }
        

        
        int thread_ID = G4Threading::G4GetThreadId();
        if (G4Threading::IsMasterThread() == true) {
            thread_ID = 0;
            //std::cout << "ActionInitialization::Build(), new thread_ID = " << thread_ID << std::endl;
        }
        
        
        //get the individual objects to be used in the threads
        PerEventObjectsManager* per_Event_Objects_Mgr = &(perEvtObjManagerVector.at(thread_ID));
        PhysicsProcessInfo* physics_Process_Info = &(physicsProcessInfoVector.at(thread_ID));
        SensitiveDetectorHelper* sensitive_Detector_Helper = &(sensitiveDetectorHelperVector.at(thread_ID));
        
        //std::cout << "Address of physics_Process_Info for thread " << G4Threading::G4GetThreadId() << "is " << physics_Process_Info << std::endl;
        //std::cout << "Address of per_Event_Objects_Mgr for thread " << G4Threading::G4GetThreadId() << "is " << per_Event_Objects_Mgr << std::endl;
        //std::cout << "Address of sensitive_Detector_Helper for thread " << G4Threading::G4GetThreadId() << "is " << sensitive_Detector_Helper << std::endl;
 
        
        PrimaryGeneratorAction* genAction = new PrimaryGeneratorAction(pset_, _genEventBroker, per_Event_Objects_Mgr);
        SetUserAction(genAction);
    
        
        Mu2eG4SteppingAction* steppingAction = new Mu2eG4SteppingAction(pset_, timeVDtimes_, _steppingCuts, _commonCuts, trajectoryControl_, mu2elimits_);
        SetUserAction(steppingAction);
    
    
        SetUserAction( new Mu2eG4StackingAction(pset_, _stackingCuts, _commonCuts) );
    
    
        TrackingAction* trackingAction = new TrackingAction(pset_, steppingAction, trajectoryControl_, mu2elimits_);
        SetUserAction(trackingAction);
    
    
        SetUserAction( new Mu2eG4RunAction(use_G4MT_, originInWorld, _physVolHelper,
                                           physics_Process_Info, trackingAction, steppingAction, sensitive_Detector_Helper) );

        
        SetUserAction( new Mu2eG4EventAction(pset_, trackingAction, steppingAction, _extMonFNALPixelSD,
                                             sensitive_Detector_Helper, _stackingCuts, _steppingCuts,
                                             _commonCuts, _genEventBroker, per_Event_Objects_Mgr,
                                             physics_Process_Info ) );
        
            }
    
    
G4VSteppingVerbose* ActionInitialization::InitializeSteppingVerbose() const
     {
         return new SteppingVerbose;
     }
    
    
    
} // end namespace mu2e
