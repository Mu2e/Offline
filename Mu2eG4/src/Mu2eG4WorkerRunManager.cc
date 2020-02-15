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
#include "G4UserWorkerThreadInitialization.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ParallelWorldProcessStore.hh"
#include "G4ParticleHPManager.hh"
#include "G4HadronicProcessStore.hh"

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
  Mu2eG4WorkerRunManager::Mu2eG4WorkerRunManager(const Mu2eG4Config::Top& conf, std::thread::id worker_ID):
    G4WorkerRunManager(),
    conf_(conf),
    m_managerInitialized(false),
    m_steppingVerbose(true),
    m_mtDebugOutput(conf.debug().mtDebugOutput()),
    rmvlevel_(conf.debug().diagLevel()),
    perThreadObjects_(std::make_unique<Mu2eG4PerThreadStorage>()),
    masterRM(nullptr),
    workerID_(worker_ID),
    mu2elimits_(conf.ResourceLimits()),
    trajectoryControl_(conf.TrajectoryControl()),
    multiStagePars_(conf),

    physicsProcessInfo_(),
    sensitiveDetectorHelper_(conf.SDConfig()),
    extMonFNALPixelSD_(),
    stackingCuts_(createMu2eG4Cuts(conf.Mu2eG4StackingOnlyCut.get<fhicl::ParameterSet>(), mu2elimits_)),
    steppingCuts_(createMu2eG4Cuts(conf.Mu2eG4SteppingOnlyCut.get<fhicl::ParameterSet>(), mu2elimits_)),
    commonCuts_(createMu2eG4Cuts(conf.Mu2eG4CommonCut.get<fhicl::ParameterSet>(), mu2elimits_))
  {
    if (m_mtDebugOutput) {
      G4cout << "WorkerRM on thread " << workerID_ << " is being created\n!";
    }

  }

  // Destructor of base is called automatically.  No need to do anything.
  Mu2eG4WorkerRunManager::~Mu2eG4WorkerRunManager(){
    if (m_mtDebugOutput) {
      G4cout << "WorkerRM on thread " << workerID_ << " is being destroyed\n!";
    }
  }


  void Mu2eG4WorkerRunManager::initializeThread(Mu2eG4MTRunManager* mRM, const G4ThreeVector& origin_in_world){

    masterRM = mRM;

    if (m_mtDebugOutput) {
      G4cout << "starting WorkerRM::initializeThread on thread: " << workerID_ << G4endl;
    }

    G4Threading::G4SetThreadId(getThreadIndex());

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
      //if(G4VSteppingVerbose::GetInstance()){std::cout << "WE HAVE STEPV2" << std::endl;}
      //WE CANNOT INSTANTIATE THIS ONE RIGHT NOW SINCE WE ALREADY HAVE ONE
      //perThreadObjects_->steppingVerbose = new SteppingVerbose();
      //SteppingVerbose* sv = perThreadObjects_->steppingVerbose;
      //if (sv)
      //    SteppingVerbose::SetInstance(sv);
    }

    // Initialize worker part of shared resources (geometry, physics)
    G4WorkerThread::BuildGeometryAndPhysicsVector();

    // Set the geometry for the worker, share from master
    G4VPhysicalVolume* worldPV = G4MTRunManager::GetMasterRunManagerKernel()->GetCurrentWorld();
    kernel->WorkerDefineWorldVolume(worldPV);

    G4TransportationManager::GetTransportationManager()->SetWorldForTracking(worldPV);

    const G4VUserDetectorConstruction* detector = masterRM->GetUserDetectorConstruction();
    G4RunManager::SetUserInitialization( const_cast<G4VUserDetectorConstruction*>(detector) );
    const_cast<G4VUserDetectorConstruction*>(detector)->ConstructSDandField();

    // Set the physics list for the worker, share from master
    physicsList = const_cast<G4VUserPhysicsList*>(masterRM->GetUserPhysicsList());
    SetUserInitialization(physicsList);

    physicsList->SetVerboseLevel(rmvlevel_);
    SetVerboseLevel(rmvlevel_);
    G4ParticleHPManager::GetInstance()->SetVerboseLevel(rmvlevel_);
    G4HadronicProcessStore::Instance()->SetVerbose(rmvlevel_);

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

    G4StateManager* stateManager = G4StateManager::GetStateManager();
    G4String currentState =  stateManager->GetStateString(stateManager->GetCurrentState());

    //we have to do this so that the state is correct for RunInitialization
    G4StateManager::GetStateManager()->SetNewState(G4State_Idle);
    if (m_mtDebugOutput) {
      G4cout << "completed WorkerRM::initializeThread on thread " << workerID_ << G4endl;
    }
  }


  void Mu2eG4WorkerRunManager::initializeUserActions(const G4ThreeVector& origin_in_world){

    userPrimaryGeneratorAction = new PrimaryGeneratorAction(conf_.debug(), perThreadObjects_.get());
    SetUserAction(userPrimaryGeneratorAction);

    steppingAction_ = new Mu2eG4SteppingAction(conf_.debug(),
                                               conf_.SDConfig().TimeVD().times(),
                                               *steppingCuts_.get(),
                                               *commonCuts_.get(),
                                               trajectoryControl_,
                                               mu2elimits_);
    SetUserAction(steppingAction_);

    SetUserAction( new Mu2eG4StackingAction(*stackingCuts_.get(),
                                            *commonCuts_.get()) );

    trackingAction_ = new TrackingAction(conf_,
                                         steppingAction_,
                                         multiStagePars_.simParticleNumberOffset(),
                                         trajectoryControl_,
                                         mu2elimits_);
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
                                         *stackingCuts_.get(),
                                         *steppingCuts_.get(),
                                         *commonCuts_.get(),
                                         perThreadObjects_.get(),
                                         &physicsProcessInfo_,
                                         origin_in_world) );

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

    m_managerInitialized = true;
  }


  void Mu2eG4WorkerRunManager::processEvent(art::Event* event){

    numberOfEventToBeProcessed = 1;
    numberOfEventProcessed = 0;
    ConstructScoringWorlds();

    eventLoopOnGoing = true;
    while(seedsQueue.size()>0)
      { seedsQueue.pop(); }
    // for each run, worker should receive at least one set of random number seeds.
    runIsSeeded = false;

    eventLoopOnGoing = true;

    ProcessOneEvent(event->id().event());

  }


} // end namespace mu2e
