//
// Override the G4RunManager class so that the Mu2e framework can drive
// the event loop.
//
// Original author Lisa Goodenough
//
//
// Notes:
//
// Implementation file for Mu2eG4WorkerRunManager

//Framework includes
#include "cetlib_except/exception.h"

//Mu2e includes
#include "Mu2eG4/inc/Mu2eG4WorkerRunManager.hh"
#include "Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Mu2eG4/inc/Mu2eG4MTRunManager.hh"
#include "Mu2eG4/inc/Mu2eG4WorkerInitialization.hh"
//#include "Mu2eG4/inc/ActionInitialization.hh"
#include "Mu2eG4/inc/SteppingVerbose.hh"
#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/Mu2eG4SteppingAction.hh"
#include "Mu2eG4/inc/Mu2eG4StackingAction.hh"
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/Mu2eG4RunAction.hh"
#include "Mu2eG4/inc/Mu2eG4EventAction.hh"
#include "Mu2eG4/inc/Mu2eG4MultiStageParameters.hh"
#include "Mu2eG4/inc/ExtMonFNALPixelSD.hh"

//C++ includes
#include <atomic>

//G4 includes
#include "G4WorkerThread.hh"
#include "G4StateManager.hh"
//#include "G4VSteppingVerbose.hh"
#include "G4UserWorkerThreadInitialization.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4VUserPhysicsList.hh"

//Other includes
#include "CLHEP/Random/JamesRandom.h"
#include <tbb/atomic.h>


using namespace std;

namespace {
    tbb::atomic<int> thread_counter{0};
    
    int get_new_thread_index() { return thread_counter++; }
    
    thread_local int s_thread_index = get_new_thread_index();
    
    int getThreadIndex() { return s_thread_index; }
    
    
}

namespace mu2e {
    
//thread_local Mu2eG4PerThreadStorage* Mu2eG4WorkerRunManager::PerThreadObjects_{nullptr};
    
    
// If the c'tor is called a second time, the c'tor of base will
// generate an exception.
Mu2eG4WorkerRunManager::Mu2eG4WorkerRunManager(const fhicl::ParameterSet& pset):
    G4WorkerRunManager(),
    pset_(pset),
    m_managerInitialized(false),
    m_userWorkerInit(true),
    m_steppingVerbose(true),
    perThreadObjects_(nullptr),
    masterRM(nullptr),
    threadID_(-10),
    mu2elimits_(pset.get<fhicl::ParameterSet>("ResourceLimits")),
    trajectoryControl_(pset.get<fhicl::ParameterSet>("TrajectoryControl")),
    multiStagePars_(pset.get<fhicl::ParameterSet>("MultiStageParameters")),
    
    physicsProcessInfo_(),
    sensitiveDetectorHelper_(pset.get<fhicl::ParameterSet>("SDConfig", fhicl::ParameterSet())),
    extMonFNALPixelSD_(),
    stackingCuts_(createMu2eG4Cuts(pset_.get<fhicl::ParameterSet>("Mu2eG4StackingOnlyCut", fhicl::ParameterSet()), mu2elimits_)),
    steppingCuts_(createMu2eG4Cuts(pset_.get<fhicl::ParameterSet>("Mu2eG4SteppingOnlyCut",fhicl::ParameterSet()), mu2elimits_)),
    commonCuts_(createMu2eG4Cuts(pset_.get<fhicl::ParameterSet>("Mu2eG4CommonCut", fhicl::ParameterSet()), mu2elimits_))
{
    SetVerboseLevel(2);
    
    
}
  
// Destructor of base is called automatically.  No need to do anything.
Mu2eG4WorkerRunManager::~Mu2eG4WorkerRunManager(){
    std::cout << "This WorkerRM is being destroyed!" << std::endl;
}
    
    
void Mu2eG4WorkerRunManager::initializePTS(Mu2eG4PerThreadStorage* pls){

    perThreadObjects_ = pls;
    threadID_ = getThreadIndex();
        
    std::cerr << "WorkerRM PLS " << threadID_ << " is Initialized!!\n";
        
}
    

void Mu2eG4WorkerRunManager::initializeThread(Mu2eG4MTRunManager* mRM, const G4ThreeVector& origin_in_world){
        
    masterRM = mRM;
    
    std::cout << "starting WorkerRM::initializeThread on thread: " << threadID_ << std::endl;
    std::cout << "test_INT = " << perThreadObjects_->test_INT << std::endl;
        
    G4Threading::G4SetThreadId(getThreadIndex());
        
    //should this be TL?
    const CLHEP::HepRandomEngine* masterEngine = masterRM->getMasterRandomEngine();
    masterRM->GetUserWorkerThreadInitialization()->SetupRNGEngine(masterEngine);
    
    //if(masterRM->GetUserWorkerInitialization())
    if(m_userWorkerInit)
    {
        //we don't have one of these
        //perThreadObjects_->workerInit = new Mu2eG4WorkerInitialization(pset_);
        //perThreadObjects_->workerInit->WorkerInitialize();
    }

    //perThreadObjects_->UserActionInit->InitializeSteppingVerbose()
    if(m_steppingVerbose) {
        if(perThreadObjects_->steppingVerbose){std::cout << "WE HAVE STEPV1" << std::endl;}
        if(G4VSteppingVerbose::GetInstance()){std::cout << "WE HAVE STEPV2" << std::endl;}
            //WE CANNOT INSTANTIATE THIS ONE RIGHT NOW SINCE WE ALREADY HAVE ONE
            //perThreadObjects_->steppingVerbose = new SteppingVerbose();
            //SteppingVerbose* sv = perThreadObjects_->steppingVerbose;
            //if (sv)
            //    SteppingVerbose::SetInstance(sv);
    }
    
    // Initialize worker part of shared resources (geometry, physics)
    G4WorkerThread::BuildGeometryAndPhysicsVector();

    // Create unique_ptr to worker run manager
    //perThreadObjects_->kernel.reset(G4WorkerRunManagerKernel::GetRunManagerKernel());
    //if (!perThreadObjects_->kernel) {
    //    perThreadObjects_->kernel.reset(new G4WorkerRunManagerKernel());
    //}
    perThreadObjects_->kernel = G4WorkerRunManagerKernel::GetRunManagerKernel();
    if (!perThreadObjects_->kernel) {
        std::cout << "Making a NEW kernel" << std::endl;
        perThreadObjects_->kernel = new G4WorkerRunManagerKernel();
    }
  
        
    // Set the geometry for the worker, share from master
    //This stuff happens in this call: localRM->Initialize();
 
    //G4RunManagerKernel* masterKernel = G4MTRunManager::GetMasterRunManagerKernel();
    //G4VPhysicalVolume* worldPV = masterKernel->GetCurrentWorld();
    G4VPhysicalVolume* worldPV = G4MTRunManager::GetMasterRunManagerKernel()->GetCurrentWorld();
    perThreadObjects_->kernel->WorkerDefineWorldVolume(worldPV);
        
    //G4TransportationManager* tM = G4TransportationManager::GetTransportationManager();
    //tM->SetWorldForTracking(worldPV);
    G4TransportationManager::GetTransportationManager()->SetWorldForTracking(worldPV);
    
    //try this once I get the ptr made in perThreadStorage
    //perThreadObjects_->userDetector = masterRM->GetUserDetectorConstruction();
    //this->G4RunManager::SetUserInitialization( perThreadObjects_->userDetector );
    //perThreadObjects_->userDetector->ConstructSDandField();
        
    const G4VUserDetectorConstruction* detector = masterRM->GetUserDetectorConstruction();
    G4RunManager::SetUserInitialization( const_cast<G4VUserDetectorConstruction*>(detector) );
    const_cast<G4VUserDetectorConstruction*>(detector)->ConstructSDandField();
    
    // Set the physics list for the worker, share from master
    std::cout << "WorkerRunManager initializing the PhysicsList" << std::endl;
    physicsList = const_cast<G4VUserPhysicsList*>(masterRM->GetUserPhysicsList());
    SetUserInitialization(physicsList);
    //const G4VUserPhysicsList* workerPhysicsList = masterRM->GetUserPhysicsList();
    //SetUserInitialization(const_cast<G4VUserPhysicsList*>(workerPhysicsList));
    
    //these two calls made in SetUserInit when called by workerRM
    //const_cast<G4VUserPhysicsList*>(perThreadObjects_->physicsList)->InitializeWorker();
    //perThreadObjects_->kernel->SetPhysics( const_cast<G4VUserPhysicsList*>(perThreadObjects_->physicsList) );
    
    //these called in G4RunManager::InitializePhysics()
    G4StateManager::GetStateManager()->SetNewState(G4State_Init);
    perThreadObjects_->kernel->InitializePhysics();
 
    const bool kernelInit = perThreadObjects_->kernel->RunInitialization();
    if (!kernelInit) {
        throw cet::exception("WorkerRUNMANAGER")
        << "Error: WorkerRunManager Geant4 kernel initialization failed!\n";
    }
    
    initializeUserActions(origin_in_world);
        
    if(masterRM->GetUserWorkerInitialization())
        { masterRM->GetUserWorkerInitialization()->WorkerStart(); }
        
    G4StateManager* stateManager = G4StateManager::GetStateManager();
    G4String currentState =  stateManager->GetStateString(stateManager->GetCurrentState());
    std::cout << "Current G4State is : " << currentState << std::endl;
        
    //Initialize();//This is really just a check at this point
        
    //we have to do this so that the state is correct for RunInitialization
    G4StateManager::GetStateManager()->SetNewState(G4State_Idle);
    std::cout << "completed WorkerRM::initializeThread on thread " << threadID_ << std::endl;
        
}
    
   
void Mu2eG4WorkerRunManager::initializeUserActions(const G4ThreeVector& origin_in_world){
        
    std::cout << "We are in WorkerRM::InitializeUserActions on thread " << threadID_ << std::endl;
    
    genAction_ = new PrimaryGeneratorAction(pset_, perThreadObjects_);
    SetUserAction(genAction_);
    
    steppingAction_ = new Mu2eG4SteppingAction(pset_,
                                               pset_.get<std::vector<double> >("SDConfig.TimeVD.times"),
                                               *steppingCuts_.get(),
                                               *commonCuts_.get(),
                                               trajectoryControl_,
                                               mu2elimits_);
    SetUserAction(steppingAction_);

    SetUserAction( new Mu2eG4StackingAction(pset_,
                                            *stackingCuts_.get(),
                                            *commonCuts_.get()) );
    
    trackingAction_ = new TrackingAction(pset_,
                                         steppingAction_,
                                         multiStagePars_.simParticleNumberOffset(),
                                         trajectoryControl_,
                                         mu2elimits_);
    SetUserAction(trackingAction_);

    SetUserAction( new Mu2eG4RunAction(pset_,
                                       origin_in_world,
                                       masterRM->getPhysVolumeHelper(),
                                       &physicsProcessInfo_,
                                       trackingAction_,
                                       steppingAction_,
                                       &sensitiveDetectorHelper_,
                                       extMonFNALPixelSD_) );
    
    SetUserAction( new Mu2eG4EventAction(pset_,
                                         trackingAction_,
                                         steppingAction_,
                                         extMonFNALPixelSD_,
                                         &sensitiveDetectorHelper_,
                                         *stackingCuts_.get(),
                                         *steppingCuts_.get(),
                                         *commonCuts_.get(),
                                         perThreadObjects_,
                                         &physicsProcessInfo_,
                                         origin_in_world) );

        /*
         if(masterRM->GetUserActionInitialization())
         { masterRM->GetNonConstUserActionInitialization()->Build(); }
         if(masterRM->GetUserWorkerInitialization())
         { masterRM->GetUserWorkerInitialization()->WorkerStart(); }
         */
}

void Mu2eG4WorkerRunManager::initializeRun(art::Event* art_event){
    
    //NOTE: we may need to add some more functionality here
    
    if (art_event->id().run() != perThreadObjects_->currentRunNumber) {
        if (perThreadObjects_->currentRunNumber != 0 && !perThreadObjects_->runTerminated) {
            //terminateRun();
            throw cet::exception("WorkerRUNMANAGER") << "Error: There is a problem with Run Numbering\n";
        }
        
        perThreadObjects_->currentRun = new G4Run();
        G4StateManager::GetStateManager()->SetNewState(G4State_GeomClosed);
        
        if(userRunAction){
            userRunAction->BeginOfRunAction(currentRun);
        }
        //if (m_pls->userRunAction) {
        //    m_pls->userRunAction->BeginOfRunAction(perThreadObjects_->currentRun);
        //}
            
        perThreadObjects_->currentRunNumber = art_event->id().run();
    }
    
        m_managerInitialized = true;
}

    
void Mu2eG4WorkerRunManager::processEvent(art::Event* event){
    
    numberOfEventToBeProcessed = 1;
    numberOfEventProcessed = 0;
    ConstructScoringWorlds();
    
    std::cout << "WorkerRM::ProcessEvent:" << event->id().event() << " on thread " << threadID_ << std::endl;
    //RunInitialization();DONE ABOVE
    
    ProcessOneEvent(event->id().event());
    TerminateOneEvent();
    
    //DoEventLoop(1);
    //RunTermination();SKIP
    
    
    //_runManager->InitializeEventLoop(num_events,macroFile,n_select);
    //_timer->Start();
    //_runManager->ProcessOneEvent(eventNumber);
    //_timer->Stop();
    
    // Accumulate time spent in G4 for all events in this run.
    //_realElapsed   += _timer->GetRealElapsed();
    //_systemElapsed += _timer->GetSystemElapsed();
    //_userElapsed   += _timer->GetUserElapsed();
    //_runManager->TerminateOneEvent();
        
    
        
}
    
    
} // end namespace mu2e
