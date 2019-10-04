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
#include "G4Threading.hh"

//art includes
#include "fhiclcpp/ParameterSet.h"



namespace mu2e {
    
    ActionInitialization::ActionInitialization(const fhicl::ParameterSet& pset,
                                               std::vector< SensitiveDetectorHelper> &sensitive_detectorhelper_vector,
                                               GenEventBroker* gen_eventbroker,
                                               PhysicalVolumeHelper* phys_volume_helper,
                                               const bool using_MT,
                                               const int num_threads,
                                               CLHEP::Hep3Vector const& origin_in_world,
                                               Mu2eG4ResourceLimits const& mu2e_limits,
                                               unsigned stage_offset_for_tracking_action
                                               )
        :
        G4VUserActionInitialization(),
        pset_(pset),
        trajectoryControl_(pset.get<fhicl::ParameterSet>("TrajectoryControl")),
        simParticlePrinter_(pset.get<fhicl::ParameterSet>("SimParticlePrinter", SimParticleCollectionPrinter::defaultPSet())),
        timeVDtimes_(pset.get<std::vector<double> >("SDConfig.TimeVD.times")),
        mu2elimits_(pset.get<fhicl::ParameterSet>("ResourceLimits")),
        _genEventBroker(gen_eventbroker),
        _physVolHelper(phys_volume_helper),
        processInfo(),
        use_G4MT_(using_MT),
        numthreads(num_threads),
        originInWorld(origin_in_world),
        mu2eLimits(mu2e_limits),
        stageOffset(stage_offset_for_tracking_action),
        sensitiveDetectorHelperVector(sensitive_detectorhelper_vector)
        {
            perEvtObjManagerVector.reserve(numthreads);
            physicsProcessInfoVector.reserve(numthreads);
            stackingCutsVector.reserve(numthreads);
            steppingCutsVector.reserve(numthreads);
            commonCutsVector.reserve(numthreads);
            
            for (int i = 0; i < numthreads ; i++) {
                perEvtObjManagerVector.emplace_back(pset_,numthreads);
                physicsProcessInfoVector.emplace_back();
                
                stackingCutsVector.emplace_back(createMu2eG4Cuts(pset_.get<fhicl::ParameterSet>("Mu2eG4StackingOnlyCut", fhicl::ParameterSet()), mu2eLimits));
                steppingCutsVector.emplace_back(createMu2eG4Cuts(pset_.get<fhicl::ParameterSet>("Mu2eG4SteppingOnlyCut", fhicl::ParameterSet()), mu2eLimits));
                commonCutsVector.emplace_back(createMu2eG4Cuts(pset_.get<fhicl::ParameterSet>("Mu2eG4CommonCut", fhicl::ParameterSet()), mu2eLimits));
                
            }

        }
    
    ActionInitialization::~ActionInitialization()
    {}

    
// used for defining only the UserRunAction for the master thread
    void ActionInitialization::BuildForMaster() const
    {
      SetUserAction( new Mu2eG4MasterRunAction(pset_, _physVolHelper, &physicsProcessInfoVector) );
    }
    
    
// used for defining user action classes for worker threads as well as for the sequential mode.
    void ActionInitialization::Build() const
    {        
    
        int thread_ID = G4Threading::G4GetThreadId();
        if (G4Threading::IsMasterThread() == true) {
            thread_ID = 0;
        }
        
        //get the individual objects to be used in the threads
        PerEventObjectsManager* per_Event_Objects_Mgr = &(perEvtObjManagerVector[thread_ID]);
        PhysicsProcessInfo* physics_Process_Info = &(physicsProcessInfoVector[thread_ID]);
        SensitiveDetectorHelper* sensitive_Detector_Helper = &(sensitiveDetectorHelperVector[thread_ID]);
        
        IMu2eG4Cut& stacking_Cuts = *stackingCutsVector[thread_ID].get();
        IMu2eG4Cut& stepping_Cuts = *steppingCutsVector[thread_ID].get();
        IMu2eG4Cut& common_Cuts = *commonCutsVector[thread_ID].get();
        
        PrimaryGeneratorAction* genAction = new PrimaryGeneratorAction(pset_, _genEventBroker, per_Event_Objects_Mgr);
        SetUserAction(genAction);
    
        Mu2eG4SteppingAction* steppingAction = new Mu2eG4SteppingAction(pset_, timeVDtimes_, stepping_Cuts, common_Cuts, trajectoryControl_, mu2elimits_);
        SetUserAction(steppingAction);
        
        SetUserAction( new Mu2eG4StackingAction(pset_, stacking_Cuts, common_Cuts) );
    
        TrackingAction* trackingAction = new TrackingAction(pset_, steppingAction, stageOffset, trajectoryControl_, mu2elimits_);
        SetUserAction(trackingAction);
        
        SetUserAction( new Mu2eG4RunAction(pset_, use_G4MT_, originInWorld, _physVolHelper,
                                           physics_Process_Info, trackingAction, steppingAction,
                                           sensitive_Detector_Helper) );

        
        SetUserAction( new Mu2eG4EventAction(pset_, trackingAction, steppingAction,
                                             sensitive_Detector_Helper, stacking_Cuts, stepping_Cuts,
                                             common_Cuts, _genEventBroker, per_Event_Objects_Mgr,
                                             physics_Process_Info,
                                             originInWorld) );
        
    }//Build()
    
    
G4VSteppingVerbose* ActionInitialization::InitializeSteppingVerbose() const
     {
         return new SteppingVerbose;
     }
    
    
} // end namespace mu2e
