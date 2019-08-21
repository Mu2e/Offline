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

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "GeometryService/inc/WorldG4.hh"

#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"

#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/preG4InitializeTasks.hh"

#include "Mu2eG4/inc/ActionInitialization.hh"
//#include "Mu2eG4/inc/Mu2eG4WorkerInitialization.hh"
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
#include "G4MTHepRandom.hh"


using namespace std;

namespace mu2e {
  
    
Mu2eG4MTRunManager* Mu2eG4MTRunManager::fMu2eMasterRM = 0;
    
// If the c'tor is called a second time, the c'tor of base will
// generate an exception.
Mu2eG4MTRunManager::Mu2eG4MTRunManager(const fhicl::ParameterSet& pset):
    G4MTRunManager(),
    pset_(pset),
    m_managerInitialized(false),
    m_runTerminated(false),
    physVolHelper_(nullptr),
    sensitiveDetectorHelper_(pset.get<fhicl::ParameterSet>("SDConfig", fhicl::ParameterSet())),
    masterRunAction_(nullptr),
    physicsList_(nullptr),
    rmvlevel_(pset.get<int>("debug.diagLevel",0))
    {
        if ( fMu2eMasterRM )
        {
            throw cet::exception("MTRUNMANAGER")
            << "Error: you are trying to create an MTRunManager when one already exists!\n";
        }
        fMu2eMasterRM = this;
    }
  
// Destructor of base is called automatically.  No need to do anything.
Mu2eG4MTRunManager::~Mu2eG4MTRunManager()
    {}
    
    
Mu2eG4MTRunManager* Mu2eG4MTRunManager::GetMasterRunManager()
    {
        return fMu2eMasterRM;
    }
    
    
void Mu2eG4MTRunManager::initializeG4(int art_runnumber)
    {
        if (m_managerInitialized) {
            std::cout << "Mu2eG4MTRunManager::initializeG4 was already done - exit";
            return;
        }
        
        SetNumberOfThreads(1);
            
        declarePhysicsAndGeometry();
            
        //Below here is similar to _runManager->RunInitialization();
        initializeKernelAndRM();
        initializeMasterRunAction();
            
        G4StateManager* stateManager = G4StateManager::GetStateManager();
        G4String currentState = stateManager->GetStateString(stateManager->GetCurrentState());
        std::cout << "Current G4State is : " << currentState << std::endl;
            
        if (GetCurrentRun()) {
            delete currentRun;
        }
        else {
            currentRun = new G4Run();
            currentRun->SetRunID(art_runnumber);
        }
            
        std::cout << "Art Run Number is: " << art_runnumber << std::endl;
        std::cout << "Current Run is " << GetCurrentRun()->GetRunID() << std::endl;
            
        //We don't need to do this here.  RunActions are owned by the worker threads
        //m_userRunAction->BeginOfRunAction(m_currentRun);
        //if(userRunAction) userRunAction->BeginOfRunAction(currentRun);
            
        currentRun->SetDCtable(DCtable);
        G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
        if(fSDM)
            { currentRun->SetHCtable(fSDM->GetHCtable()); }
        std::ostringstream oss;
        //WHAT LIBRARY TO USE?
        //G4Random::saveFullState(oss);
        //randomNumberStatusForThisRun = oss.str();
        //currentRun->SetRandomNumberStatus(randomNumberStatusForThisRun);
        
            
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
/////////////////////////////////////////////////////////////////////////////////////
        //RANDOM NUMBER SEEDING: NEED TO FIGURE THIS OUT
        //taken from G4MTRunManager::InitializeEventLoop()
        int numevents = 1000;
        int numworkers = 1;
        SetNumberOfEventsToBeProcessed(numevents);
        //numberOfEventToBeProcessed = 1000;//this needs to be filled from the # of events we are processing
        numberOfEventProcessed = 0;
        
        nSeedsUsed = 0;
        nSeedsFilled = 0;
        
        if(verboseLevel>0)
        { timer->Start(); }
   
        //n_select_msg = -1;
        //selectMacro = "";

        //initialize seeds
        //If user did not implement InitializeSeeds,
        // use default: nSeedsPerEvent seeds per event
        //eventModuloDef = 0 by default
        if( eventModuloDef > 0 ) {
            eventModulo = eventModuloDef;
            if(eventModulo > numberOfEventToBeProcessed/numworkers) {
                eventModulo = numberOfEventToBeProcessed/numworkers;
                if(eventModulo<1) eventModulo =1;
                
                G4ExceptionDescription msgd;
                msgd << "Event modulo is reduced to " << eventModulo
                     << " to distribute events to all threads.";
                G4Exception("G4MTRunManager::InitializeEventLoop()", "Run10035", JustWarning, msgd);
            }
        }
        else {
            eventModulo = int(std::sqrt(double(numberOfEventToBeProcessed/numworkers)));
            if(eventModulo<1) eventModulo = 1;
        }
        
        if ( InitializeSeeds(numevents) == false && numevents>0 ) {
            
            G4RNGHelper* helper = G4RNGHelper::GetInstance();
            switch(seedOncePerCommunication)
            {
                case 0://default value
                    nSeedsFilled = numevents;
                    break;
                case 1:
                    nSeedsFilled = numworkers;
                    break;
                case 2:
                    nSeedsFilled = numevents/eventModulo + 1;
                    break;
                default:
                    G4ExceptionDescription msgd;
                    msgd << "Parameter value <" << seedOncePerCommunication
                         << "> of seedOncePerCommunication is invalid. It is reset to 0." ;
                    G4Exception("G4MTRunManager::InitializeEventLoop()", "Run10036", JustWarning, msgd);
                    seedOncePerCommunication = 0;
                    nSeedsFilled = numevents;
            }

            // Generates up to nSeedsMax seed pairs only.
            if(nSeedsFilled>nSeedsMax) nSeedsFilled=nSeedsMax;
            const_cast<CLHEP::HepRandomEngine*>(getMasterRandomEngine())->flatArray(nSeedsPerEvent*nSeedsFilled,randDbl);
            helper->Fill(randDbl,nSeedsFilled,numevents,nSeedsPerEvent);
        }

      

    }
    
    
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
            allMu2e = (new WorldMaker<Mu2eWorld>(make_unique<Mu2eWorld>(pset_, &sensitiveDetectorHelper_),
                                                 make_unique<ConstructMaterials>(pset_)));
        }
        else {
            throw cet::exception("CONFIG")
            << "Error: You are trying to run in MT mode without the Standard Mu2e Detector!\n";
        }
        
        //G4ThreeVector originInWorld = GeomHandle<WorldG4>()->mu2eOriginInWorld();
        
        preG4InitializeTasks(pset_);
        
        physicsList_ = physicsListDecider(pset_);
        
        physicsList_->SetVerboseLevel(rmvlevel_);
        SetVerboseLevel(rmvlevel_);
        G4ParticleHPManager::GetInstance()->SetVerboseLevel(rmvlevel_);
        G4HadronicProcessStore::Instance()->SetVerbose(rmvlevel_);
        
        SetUserInitialization(allMu2e);
        SetUserInitialization(physicsList_);
        //SetUserInitialization( new Mu2eG4WorkerInitialization(pset_) );//This is NEW code, not sure we need it.
        
        
        /*ActionInitialization* actionInit = new ActionInitialization(pset_,
                                                                    sensitiveDetectorHelper_,
                                                                    &_genEventBroker,
                                                                    &physVolHelper_,
                                                                    originInWorld
                                                                    );
        */
        //in MT mode, this is where BuildForMaster is called for master thread
        // in sequential mode, this is where Build() is called for main thread
        //        SetUserInitialization(actionInit);
        
    }
  
    
void Mu2eG4MTRunManager::initializeMasterRunAction()
    {
        masterRunAction_ = new Mu2eG4MasterRunAction(pset_, physVolHelper_);
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
    
        /*if (m_userRunAction) {
            m_userRunAction->EndOfRunAction(m_currentRun);
            delete m_userRunAction;
            m_userRunAction = nullptr;
        }*/
        if ((G4MTRunManager::GetMTMasterRunManagerKernel()!=nullptr) && !m_runTerminated) {
            G4MTRunManager::GetMTMasterRunManagerKernel()->RunTermination();
        }
        m_runTerminated = true;
}
    
    
void Mu2eG4MTRunManager::Test_Func(int in)
    {
        std::cout << "FROM MASTER_RM::Test_Func, the integer is " << in << std::endl;
        physicsList_->DumpList();
        
    }
    

////////////////////////// OLD STUFF BELOW HERE //////////////////////////
    

   
    //we need control of the event loop in order to correctly break up the stages
    //to fit within the art framework.  we were getting a hang using the G4MTRunManager
    //version of RunTermination due to the call of WaitForEndEventLoopWorkers()
    void Mu2eG4MTRunManager::mu2eG4RunTermination()
    {
        
        if ( verboseLevel > 0 ) {
            G4cout << __func__ << " called" << G4endl;
        }
        
        G4RunManager::TerminateEventLoop();
        G4RunManager::RunTermination();
    }
    
    
    
    
    
    
    
} // end namespace mu2e
