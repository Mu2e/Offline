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
#include "Offline/Mu2eG4/inc/Mu2eG4WorkerRunManager.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4MTRunManager.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4SteppingVerbose.hh"
#include "Offline/Mu2eG4/inc/WorldMaker.hh"
#include "Offline/Mu2eG4/inc/physicsListDecider.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4PrimaryGeneratorAction.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4SteppingAction.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4StackingAction.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4TrackingAction.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4RunAction.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4EventAction.hh"
#include "Offline/Mu2eG4/inc/ExtMonFNALPixelSD.hh"

//G4 includes
#include "Geant4/G4WorkerThread.hh"
#include "Geant4/G4StateManager.hh"
#include "Geant4/G4UserWorkerThreadInitialization.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4TransportationManager.hh"
#include "Geant4/G4VUserPhysicsList.hh"
#include "Geant4/G4ParallelWorldProcessStore.hh"
#if G4VERSION>4106
#include "Geant4/G4HadronicParameters.hh"
#else
#include "Geant4/G4ParticleHPManager.hh"
#include "Geant4/G4HadronicProcessStore.hh"
#endif
//Other includes
#include <string>
#include <atomic>
#include "CLHEP/Random/JamesRandom.h"

using namespace std;

namespace {
  atomic<int> thread_counter{0};
  int get_new_thread_index() { return thread_counter++; }
  thread_local int s_thread_index = get_new_thread_index();
  int getThreadIndex() { return s_thread_index; }
}

namespace mu2e {


  // If the c'tor is called a second time, the c'tor of base will
  // generate an exception.
  Mu2eG4WorkerRunManager::Mu2eG4WorkerRunManager(const Mu2eG4Config::Top& conf, const Mu2eG4IOConfigHelper& ioconf, thread::id worker_ID):
    G4WorkerRunManager(),
    conf_(conf),
    m_managerInitialized(false),
    m_steppingVerbose(true),
    m_mtDebugOutput(conf.debug().mtDebugOutput()),
    rmvlevel_(conf.debug().diagLevel()),
    salt_(conf.salt()),
    perThreadObjects_(make_unique<Mu2eG4PerThreadStorage>(ioconf)),
    masterRM(nullptr),
    workerID_(worker_ID),
    physicsProcessInfo_(),
    sensitiveDetectorHelper_(conf.SDConfig()),
    extMonFNALPixelSD_()
  {
    if (m_mtDebugOutput > 0) {
      G4cout << "WorkerRM on thread " << workerID_ << " is being created\n!";
    }
  }

  // Destructor of base is called automatically.  No need to do anything.
  Mu2eG4WorkerRunManager::~Mu2eG4WorkerRunManager(){

    if (m_mtDebugOutput > 0) {
      G4cout << "WorkerRM on thread " << workerID_ << " is being destroyed\n!";
    }
  }


  void Mu2eG4WorkerRunManager::initializeThread(Mu2eG4MTRunManager* mRM, const G4ThreeVector& origin_in_world){

    masterRM = mRM;

    if (m_mtDebugOutput > 0) {
      G4cout << "starting WorkerRM::initializeThread on thread: " << workerID_ << G4endl;
    }

    G4Threading::G4SetThreadId(getThreadIndex());

    const CLHEP::HepRandomEngine* masterEngine = masterRM->getMasterRandomEngine();
    masterRM->GetUserWorkerThreadInitialization()->SetupRNGEngine(masterEngine);

    // Initialize worker part of shared resources (geometry, physics)
    G4WorkerThread::BuildGeometryAndPhysicsVector();

    // Set the geometry for the worker, share from master
    G4VPhysicalVolume* worldPV = G4MTRunManager::GetMasterRunManagerKernel()->GetCurrentWorld();
    kernel->WorkerDefineWorldVolume(worldPV);
    kernel->SetNumberOfParallelWorld(G4MTRunManager::GetMasterRunManagerKernel()->GetNumberOfParallelWorld());

    G4TransportationManager::GetTransportationManager()->SetWorldForTracking(worldPV);

    const G4VUserDetectorConstruction* detector = masterRM->GetUserDetectorConstruction();
    G4RunManager::SetUserInitialization( const_cast<G4VUserDetectorConstruction*>(detector) );
    const_cast<G4VUserDetectorConstruction*>(detector)->ConstructSDandField();

    // Set the physics list for the worker, share from master
    physicsList = const_cast<G4VUserPhysicsList*>(masterRM->GetUserPhysicsList());
    SetUserInitialization(physicsList);

    physicsList->SetVerboseLevel(rmvlevel_);
    SetVerboseLevel(rmvlevel_);
#if G4VERSION>4106
    G4HadronicParameters::Instance()->SetVerboseLevel(rmvlevel_);
#else
    G4ParticleHPManager::GetInstance()->SetVerboseLevel(rmvlevel_);
    G4HadronicProcessStore::Instance()->SetVerbose(rmvlevel_);
#endif
    //these called in G4RunManager::InitializePhysics()
    G4StateManager::GetStateManager()->SetNewState(G4State_Init);
    kernel->InitializePhysics();

    const bool kernelInit = kernel->RunInitialization();
    if (!kernelInit) {
      throw cet::exception("WorkerRUNMANAGER")
        << "Error: WorkerRunManager Geant4 kernel initialization failed!\n";
    }

    initializeUserActions(origin_in_world);

    G4StateManager* stateManager = G4StateManager::GetStateManager();
    G4String currentState =  stateManager->GetStateString(stateManager->GetCurrentState());

    //we have to do this so that the state is correct for RunInitialization
    G4StateManager::GetStateManager()->SetNewState(G4State_Idle);
    if (m_mtDebugOutput > 0) {
      G4cout << "completed WorkerRM::initializeThread on thread " << workerID_ << G4endl;
    }
  }


  void Mu2eG4WorkerRunManager::initializeUserActions(const G4ThreeVector& origin_in_world){

    userPrimaryGeneratorAction = new Mu2eG4PrimaryGeneratorAction(conf_.debug(), perThreadObjects_.get());
    SetUserAction(userPrimaryGeneratorAction);

    steppingAction_ = new Mu2eG4SteppingAction(conf_.debug(),
                                               conf_.SDConfig().TimeVD().times(),
                                               *perThreadObjects_->steppingCuts,
                                               *perThreadObjects_->commonCuts,
                                               perThreadObjects_->ioconf.trajectoryControl(),
                                               perThreadObjects_->ioconf.mu2elimits());
    SetUserAction(steppingAction_);

    SetUserAction( new Mu2eG4StackingAction(*perThreadObjects_->stackingCuts,
                                            *perThreadObjects_->commonCuts) );

    trackingAction_ = new Mu2eG4TrackingAction(conf_,
                                         steppingAction_,
                                         perThreadObjects_.get());
    SetUserAction(trackingAction_);

    SetUserAction( new Mu2eG4RunAction(conf_.debug(),
                                       origin_in_world,
                                       masterRM->getPhysVolumeHelper(),
                                       &physicsProcessInfo_,
                                       trackingAction_,
                                       steppingAction_,
                                       &sensitiveDetectorHelper_) );

    SetUserAction( new Mu2eG4EventAction(conf_,
                                         trackingAction_,
                                         steppingAction_,
                                         &sensitiveDetectorHelper_,
                                         perThreadObjects_.get(),
                                         &physicsProcessInfo_,
                                         origin_in_world) );

  }


  void Mu2eG4WorkerRunManager::initializeRun(art::Event* const art_event){

    perThreadObjects_->currentRunNumber = art_event->id().run();

    ConstructScoringWorlds();

    //following taken from G4WorkerRunManager::RunInitialization()
    fakeRun = false;
    if(!(kernel->RunInitialization(fakeRun))) return;

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

    ostringstream oss;
    G4Random::saveFullState(oss);
    randomNumberStatusForThisRun = oss.str();
    currentRun->SetRandomNumberStatus(randomNumberStatusForThisRun);

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
        ostringstream os;
        os << "run" << currentRun->GetRunID();
        fileN = os.str();
      }
      StoreRNGStatus(fileN);
    }

    runAborted = false;
    numberOfEventProcessed = 0;
    m_managerInitialized = true;

  }


  void Mu2eG4WorkerRunManager::processEvent(const art::EventID& evtID){

    numberOfEventToBeProcessed = 1;
    runIsSeeded = false;
    eventLoopOnGoing = true;

    //below code is from ProcessOneEvent(i_event);
    currentEvent = generateEvt(evtID);

    if(eventLoopOnGoing) {
      eventManager->ProcessOneEvent(currentEvent);
      AnalyzeEvent(currentEvent);
      UpdateScoring();
    }
  }


  G4Event* Mu2eG4WorkerRunManager::generateEvt(const art::EventID& evtID){

    G4Event* anEvent = new G4Event(evtID.event());
    G4bool eventHasToBeSeeded = true;

    eventLoopOnGoing = masterRM->SetUpEvent();
    runIsSeeded = true;

    if(!eventLoopOnGoing)
      {
        delete anEvent;
        return 0;
      }

    if(eventHasToBeSeeded)
      {
        string msg = "r" + to_string(evtID.run())
          + "s" + to_string(evtID.subRun())
          + "e" + to_string(evtID.event()) + salt_;
        std::hash<string> hf;
        long rn1 = hf(msg+"1") & 0xFFFFFFFF;
        long rn2 = hf(msg+"2") & 0xFFFFFFFF;
        long seeds[3] = { rn1, rn2, 0 };
        G4Random::setTheSeeds(seeds,-1);
        runIsSeeded = true;

        if(m_mtDebugOutput > 1) {
          G4cout << "--> Event " << anEvent->GetEventID() << " starts with initial seeds ("
                 << rn1 << "," << rn2 << ")." << G4endl;
        }

      }

    //This is the filename base constructed from run and event
    const auto filename = [&] {
      ostringstream os;
      os << "run" << currentRun->GetRunID() << "evt" << anEvent->GetEventID();
      return os.str();
    };

    G4bool RNGstatusReadFromFile = false;
    if ( readStatusFromFile ) {
      //Build full path of RNG status file for this event
      ostringstream os;
      os << filename() << ".rndm";
      const G4String& randomStatusFile = os.str();
      ifstream ifile(randomStatusFile.c_str());
      if ( ifile ) { //File valid and readable
        RNGstatusReadFromFile = true;
        G4Random::restoreEngineStatus(randomStatusFile.c_str());
      }
    }


    if(storeRandomNumberStatusToG4Event==1 || storeRandomNumberStatusToG4Event==3) {
      ostringstream oss;
      G4Random::saveFullState(oss);
      randomNumberStatusForThisEvent = oss.str();
      anEvent->SetRandomNumberStatus(randomNumberStatusForThisEvent);
    }

    if(storeRandomNumberStatus && ! RNGstatusReadFromFile ) { //If reading from file, avoid to rewrite the same
      G4String fileN = "currentEvent";
      if ( rngStatusEventsFlag ) {
        fileN = filename();
      }
      StoreRNGStatus(fileN);
    }

    userPrimaryGeneratorAction->GeneratePrimaries(anEvent);
    return anEvent;

  }

} // end namespace mu2e
