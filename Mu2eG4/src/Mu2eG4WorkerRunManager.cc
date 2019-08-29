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

//G4 includes
#include "G4WorkerThread.hh"
#include "G4StateManager.hh"
//#include "G4VSteppingVerbose.hh"
#include "G4UserWorkerThreadInitialization.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ParallelWorldProcessStore.hh"
#include "G4MTHepRandom.hh"

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
    
    
// If the c'tor is called a second time, the c'tor of base will
// generate an exception.
    Mu2eG4WorkerRunManager::Mu2eG4WorkerRunManager(const fhicl::ParameterSet& pset, std::string worker_ID):
    G4WorkerRunManager(),
    pset_(pset),
    m_managerInitialized(false),
    m_userWorkerInit(true),
    m_steppingVerbose(true),
    perThreadObjects_(std::make_unique<Mu2eG4PerThreadStorage>(pset)),
    masterRM(nullptr),
    workerID_(worker_ID),
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
    //SetVerboseLevel(2);
}
  
// Destructor of base is called automatically.  No need to do anything.
Mu2eG4WorkerRunManager::~Mu2eG4WorkerRunManager(){
    std::cout << "This WorkerRM is being destroyed!" << std::endl;
}
    
    
void Mu2eG4WorkerRunManager::initializePTS(Mu2eG4PerThreadStorage* pls){

    //perThreadObjects_ = pls;
    //workerID_ = getThreadIndex();
    //std::cerr << "WorkerRM PLS " << workerID_ << " is Initialized!!\n";
        
}
    

void Mu2eG4WorkerRunManager::initializeThread(Mu2eG4MTRunManager* mRM, const G4ThreeVector& origin_in_world){
        
    masterRM = mRM;
    
    std::cout << "starting WorkerRM::initializeThread on thread: " << workerID_ << std::endl;
    //std::cout << "test_INT = " << perThreadObjects_->test_INT << std::endl;
        
    G4Threading::G4SetThreadId(getThreadIndex());
        
    //should this be TL?
    const CLHEP::HepRandomEngine* masterEngine = masterRM->getMasterRandomEngine();
    masterRM->GetUserWorkerThreadInitialization()->SetupRNGEngine(masterEngine);
    
    //perThreadObjects_->UserActionInit->InitializeSteppingVerbose()
    if(m_steppingVerbose) {
        
        
        //if(masterRM->GetUserActionInitialization())
        //{
        //    G4VSteppingVerbose* sv = masterRM->GetUserActionInitialization()->InitializeSteppingVerbose();
        //    if ( sv ) { G4VSteppingVerbose::SetInstance(sv); }
        //}
        
        //if(perThreadObjects_->steppingVerbose){std::cout << "WE HAVE STEPV1" << std::endl;}
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
//    perThreadObjects_->kernel = G4WorkerRunManagerKernel::GetRunManagerKernel();
//    if (!perThreadObjects_->kernel) {
//        std::cout << "Making a NEW kernel" << std::endl;
//        perThreadObjects_->kernel = new G4WorkerRunManagerKernel();
//    }
  
        
    // Set the geometry for the worker, share from master
    //This stuff happens in this call: localRM->Initialize();
    //WHAT IF CALLED wrm->Initialze();?
 
    //G4RunManagerKernel* masterKernel = G4MTRunManager::GetMasterRunManagerKernel();
    //G4VPhysicalVolume* worldPV = masterKernel->GetCurrentWorld();
    G4VPhysicalVolume* worldPV = G4MTRunManager::GetMasterRunManagerKernel()->GetCurrentWorld();
    kernel->WorkerDefineWorldVolume(worldPV);
        
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
    
  
    //these called in G4RunManager::InitializePhysics()
    
    G4StateManager::GetStateManager()->SetNewState(G4State_Init);
    kernel->InitializePhysics();
    
    //WHY IS THIS DONE HERE???????????????????????????
    const bool kernelInit = kernel->RunInitialization();
    if (!kernelInit) {
        throw cet::exception("WorkerRUNMANAGER")
        << "Error: WorkerRunManager Geant4 kernel initialization failed!\n";
    }
    
    initializeUserActions(origin_in_world);
    
    //if(masterRM->GetUserWorkerInitialization())
    //    { masterRM->GetUserWorkerInitialization()->WorkerStart(); }
        
    G4StateManager* stateManager = G4StateManager::GetStateManager();
    G4String currentState =  stateManager->GetStateString(stateManager->GetCurrentState());
    std::cout << "Current G4State is : " << currentState << std::endl;
        
    //Initialize();//This is really just a check at this point
        
    //we have to do this so that the state is correct for RunInitialization
    G4StateManager::GetStateManager()->SetNewState(G4State_Idle);
    std::cout << "completed WorkerRM::initializeThread on thread " << workerID_ << std::endl;
}
    
   
void Mu2eG4WorkerRunManager::initializeUserActions(const G4ThreeVector& origin_in_world){
        
    std::cout << "We are in WorkerRM::InitializeUserActions on thread " << workerID_ << std::endl;
    
    userPrimaryGeneratorAction = new PrimaryGeneratorAction(pset_, perThreadObjects_.get());
    SetUserAction(userPrimaryGeneratorAction);
    
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
                                         perThreadObjects_.get(),
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
    
    
    if (art_event->id().run() != perThreadObjects_->currentRunNumber) {
        if (perThreadObjects_->currentRunNumber != 0 && !perThreadObjects_->runTerminated) {
            //terminateRun();
            throw cet::exception("WorkerRUNMANAGER") << "Error: There is a problem with Run Numbering\n";
        }
    
        
    //following taken from G4WorkerRunManager::RunInitialization()
    fakeRun = false;
    if(!(kernel->RunInitialization(fakeRun))) return;

    runAborted = false;
    numberOfEventProcessed = 0;

    CleanUpPreviousEvents();
    if(currentRun) delete currentRun;
    currentRun = 0;

    if(fGeometryHasBeenDestroyed) G4ParallelWorldProcessStore::GetInstance()->UpdateWorlds();
    if(userRunAction) currentRun = userRunAction->GenerateRun();

    if(!currentRun) currentRun = new G4Run();
    currentRun->SetRunID(runIDCounter);
    currentRun->SetNumberOfEventToBeProcessed(numberOfEventToBeProcessed);
    currentRun->SetDCtable(DCtable);
    
    G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
    
    if(fSDM) currentRun->SetHCtable(fSDM->GetHCtable());
    
    //AGAIN, WHAT LIBRARY?
    //std::ostringstream oss;
    //G4Random::saveFullState(oss);
    //randomNumberStatusForThisRun = oss.str();
    //currentRun->SetRandomNumberStatus(randomNumberStatusForThisRun);

    for(G4int i_prev=0;i_prev<n_perviousEventsToBeStored;i_prev++) {
        previousEvents->push_back((G4Event*)0);
    }

    if(printModulo>=0 || verboseLevel>0) {
        G4cout << "### Run " << currentRun->GetRunID() << " starts." << G4endl;
    }
    
    if(userRunAction) userRunAction->BeginOfRunAction(currentRun);
    
    if(storeRandomNumberStatus) {
        G4String fileN = "currentRun";
        if ( rngStatusEventsFlag ) {
            std::ostringstream os;
            os << "run" << currentRun->GetRunID();
            fileN = os.str();

        }
        StoreRNGStatus(fileN);
    }
    
    }
    
    perThreadObjects_->currentRunNumber = art_event->id().run();
/////////////////////////////////////////////
    
  
    
/*
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
*/
        m_managerInitialized = true;
}

    
void Mu2eG4WorkerRunManager::processEvent(art::Event* event){
    
    numberOfEventToBeProcessed = 1;
    numberOfEventProcessed = 0;
    ConstructScoringWorlds();
    
    std::cout << "WorkerRM::ProcessEvent:" << event->id().event() << " on thread " << workerID_ << std::endl;
    //RunInitialization();DONE ABOVE
    //G4WorkerRunManager::DoEventLoop(1);DOESN'T WORK
    
    eventLoopOnGoing = true;
    while(seedsQueue.size()>0)
        { seedsQueue.pop(); }
    // for each run, worker should receive at least one set of random number seeds.
    runIsSeeded = false;
    
    eventLoopOnGoing = true;
    
    ProcessOneEvent(event->id().event());
    //TerminateOneEvent();
     
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
