//
// Override the G4RunManager class so that the Mu2e framework can drive
// the event loop.
//
// Original author Lisa Goodenough
//
//
// Notes:
//
// Implementation file for Mu2eG4MTRunManager

//Framework includes
#include "cetlib_except/exception.h"

//Mu2e includes
#include "Mu2eG4/inc/Mu2eG4MTRunManager.hh"
#include "SeedService/inc/SeedService.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/preG4InitializeTasks.hh"
#include "Mu2eG4/inc/ActionInitialization.hh"
#include "Mu2eG4/inc/Mu2eG4MasterRunAction.hh"

//G4 includes
#include "G4Timer.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ParticleHPManager.hh"
#include "G4HadronicProcessStore.hh"
#include "G4StateManager.hh"
#include "G4GeometryManager.hh"
#include "G4UserWorkerThreadInitialization.hh"
#include "G4MTRunManagerKernel.hh"
#include "G4VUserPhysicsList.hh"
#include "G4SDManager.hh"

using namespace std;

namespace {G4Mutex setUpEventMutex = G4MUTEX_INITIALIZER;}

namespace mu2e {

  // If the c'tor is called a second time, the c'tor of base will
  // generate an exception.
  Mu2eG4MTRunManager::Mu2eG4MTRunManager(const Mu2eG4Config::Top& conf):
    G4MTRunManager(),
    conf_(conf),
    m_managerInitialized(false),
    m_runTerminated(false),
    physVolHelper_(nullptr),
    sensitiveDetectorHelper_(conf.SDConfig()),
    masterRunAction_(nullptr),
    physicsList_(nullptr),
    rmvlevel_(conf.debug().diagLevel())
  {
    const_cast<CLHEP::HepRandomEngine*>(getMasterRandomEngine())->setSeed(art::ServiceHandle<SeedService>()->getSeed(),0);
  }

  // Destructor of base is called automatically.  No need to do anything.
  Mu2eG4MTRunManager::~Mu2eG4MTRunManager()
  {}


  void Mu2eG4MTRunManager::initializeG4(int art_runnumber)
  {
    if (m_managerInitialized) {
      G4cout << "Mu2eG4MTRunManager::initializeG4 was already done - exit" << G4endl;
      return;
    }

    //this is the number of G4 Worker Threads.  It allows us to use the G4MT RunManager without having to create
    // extra worker threads
    SetNumberOfThreads(1);

    declarePhysicsAndGeometry();

    //Below here is similar to _runManager->RunInitialization();
    initializeKernelAndRM();
    initializeMasterRunAction();

    G4StateManager* stateManager = G4StateManager::GetStateManager();
    G4String currentState = stateManager->GetStateString(stateManager->GetCurrentState());
    G4cout << "Current G4State is : " << currentState << G4endl;

    if (GetCurrentRun()) {
      delete currentRun;
    }
    else {
      currentRun = new G4Run();
      currentRun->SetRunID(art_runnumber);
    }

    G4cout << "Art Run Number is: " << art_runnumber << G4endl;
    G4cout << "Current Run is " << GetCurrentRun()->GetRunID() << G4endl;

    currentRun->SetDCtable(DCtable);
    G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
    if(fSDM)
      { currentRun->SetHCtable(fSDM->GetHCtable()); }
    std::ostringstream oss;
    G4Random::saveFullState(oss);
    randomNumberStatusForThisRun = oss.str();
    currentRun->SetRandomNumberStatus(randomNumberStatusForThisRun);

    if(storeRandomNumberStatus) {
      G4String fileN = "currentRun";
      if ( rngStatusEventsFlag ) {
        std::ostringstream os;
        os << "run" << currentRun->GetRunID();
        fileN = os.str();
      }
      StoreRNGStatus(fileN);
    }

    /////////////////////////////////////////////////////////////////////////////////////
    //RANDOM NUMBER SEEDING
    SetEventModulo(1);//this sets eventModuloDef
    int seedbunchsize = 10000;
    int numworkers = 1;
    numberOfEventToBeProcessed = std::numeric_limits<int>::max();
    numberOfEventProcessed = 0;

    nSeedsUsed = 0;
    nSeedsFilled = 0;

    if(verboseLevel>0)
      { timer->Start(); }

    //if we are using G4's seed filling scheme
    if ( InitializeSeeds(seedbunchsize) == false && seedbunchsize>0 ) {

      G4RNGHelper* helper = G4RNGHelper::GetInstance();
      switch(seedOncePerCommunication)
        {
        case 0://default value
          nSeedsFilled = seedbunchsize;
          break;
        case 1:
          nSeedsFilled = numworkers;
          break;
        case 2:
          nSeedsFilled = seedbunchsize/eventModulo + 1;
          break;
        default:
          G4ExceptionDescription msgd;
          msgd << "Parameter value <" << seedOncePerCommunication
               << "> of seedOncePerCommunication is invalid. It is reset to 0." ;
          G4Exception("G4MTRunManager::InitializeEventLoop()", "Run10036", JustWarning, msgd);
          seedOncePerCommunication = 0;
          nSeedsFilled = seedbunchsize;
        }

      // Generates up to nSeedsMax seed pairs only.
      if(nSeedsFilled>nSeedsMax) nSeedsFilled=nSeedsMax;
      const_cast<CLHEP::HepRandomEngine*>(getMasterRandomEngine())->flatArray(nSeedsPerEvent*nSeedsFilled,randDbl);
      helper->Fill(randDbl,nSeedsFilled,seedbunchsize,nSeedsPerEvent);
    }
  }//Mu2eG4MTRunManager::initializeG4


  // this function is a protected member of G4MTRunManager but we need to access it
  // from Mu2eG4_module, so we must make it public here
  void Mu2eG4MTRunManager::initializeKernelAndRM()
  {
    G4RunManager::Initialize();
    G4MTRunManager::GetMTMasterRunManagerKernel()->SetUpDecayChannels();//note, this is usually done in
    //InitializeEventLoop

    if ( userWorkerThreadInitialization == 0 )
      { userWorkerThreadInitialization = new G4UserWorkerThreadInitialization(); }

    G4bool cond = ConfirmBeamOnCondition();
    if (cond) {

      G4MTRunManager::ConstructScoringWorlds();
      if (G4MTRunManager::GetMTMasterRunManagerKernel()->RunInitialization()) {
        m_managerInitialized = true;
      } else {
        throw cet::exception("G4MTRunManagerKernel initialization failed!");
      }
      SetRunIDCounter(0);

    } else {
      throw cet::exception("G4RunManager ConfirmBeamOnCondition failed!");
    }

  }

  void Mu2eG4MTRunManager::declarePhysicsAndGeometry()
  {
    G4VUserDetectorConstruction* allMu2e;

    if ((art::ServiceHandle<GeometryService>())->isStandardMu2eDetector()) {
      allMu2e = (new WorldMaker<Mu2eWorld>(make_unique<Mu2eWorld>(conf_, &sensitiveDetectorHelper_),
                                           make_unique<ConstructMaterials>(conf_)));
    }
    else {
      throw cet::exception("CONFIG")
        << "Error: You are trying to run in MT mode without the Standard Mu2e Detector!\n";
    }

    preG4InitializeTasks(conf_);

    physicsList_ = physicsListDecider(conf_.physics(), conf_.debug());

    physicsList_->SetVerboseLevel(rmvlevel_);
    SetVerboseLevel(rmvlevel_);
    G4ParticleHPManager::GetInstance()->SetVerboseLevel(rmvlevel_);
    G4HadronicProcessStore::Instance()->SetVerbose(rmvlevel_);

    SetUserInitialization(allMu2e);
    SetUserInitialization(physicsList_);
  }


  void Mu2eG4MTRunManager::initializeMasterRunAction()
  {
    masterRunAction_ = new Mu2eG4MasterRunAction(conf_.debug().diagLevel(), physVolHelper_);
    SetUserAction( masterRunAction_ );
    masterRunAction_->MasterBeginRunAction();
  }


  void Mu2eG4MTRunManager::stopG4()
  {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4StateManager::GetStateManager()->SetNewState(G4State_Quit);

    if (!m_runTerminated) {
      terminateRun();
    }
  }


  void Mu2eG4MTRunManager::terminateRun() {

    masterRunAction_->MasterEndRunAction();

    if ((G4MTRunManager::GetMTMasterRunManagerKernel()!=nullptr) && !m_runTerminated) {
      std::cerr << "CALLING RunTermination() from MTRunManager\n";
      G4RunManager::TerminateEventLoop();
      G4RunManager::RunTermination();
    }
    m_runTerminated = true;
  }

  G4bool Mu2eG4MTRunManager::SetUpAnEvent(G4Event* evt, long& s1, long& s2, long& s3,
                                          G4bool reseedRequired) {

    G4AutoLock l(&setUpEventMutex);

    if( numberOfEventProcessed < numberOfEventToBeProcessed ) {

      if(reseedRequired) {
        G4RNGHelper* helper = G4RNGHelper::GetInstance();
        G4int idx_rndm = nSeedsPerEvent*(evt->GetEventID()-1);
        s1 = helper->GetSeed(idx_rndm);
        s2 = helper->GetSeed(idx_rndm+1);
        if(nSeedsPerEvent==3) s3 = helper->GetSeed(idx_rndm+2);
        nSeedsUsed++;
        //G4cout << "nSeedsUsed = " << nSeedsUsed << ", nSeedsFilled = " << nSeedsFilled << "\n";
        if(nSeedsUsed==nSeedsFilled) {
          RefillSeeds();
          //G4cout << "Refilling Seeds\n";
        }
      }
      numberOfEventProcessed++;
      return true;
    }

    return false;
  }


} // end namespace mu2e
