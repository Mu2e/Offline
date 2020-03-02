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
#include "G4MTRunManagerKernel.hh"
#include "G4UserWorkerThreadInitialization.hh"

//G4 includes
#include "G4Timer.hh"

//C++ includes
#include <thread>
#include <sstream>


using namespace std;

namespace mu2e {
  
    
    Mu2eG4MTRunManager* Mu2eG4MTRunManager::fMu2eMasterRM = 0;
    
    // If the c'tor is called a second time, the c'tor of base will
    // generate an exception.
    Mu2eG4MTRunManager::Mu2eG4MTRunManager():
        G4MTRunManager()
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
    
    // this function is a protected member of G4MTRunManager but we need to access it
    // from Mu2eG4_module, so we must make it public here
    void Mu2eG4MTRunManager::Mu2eG4Initialize(G4int n_event)
    {
        G4RunManager::Initialize();
        
        //BeamOn(0);
        if(n_event<=0) { fakeRun = true; }
        else { fakeRun = false; }
        G4bool cond = ConfirmBeamOnCondition();
        if(cond)
        {
            numberOfEventToBeProcessed = n_event;
            numberOfEventProcessed = 0;
            G4MTRunManager::ConstructScoringWorlds();
            RunInitialization();
            
            //Instead of calling DoEventLoop(n_event);
            //we call the only piece that happens for 0 events
            //and an MT RunManager
            InitializeEventLoop(n_event);
            
            RunTermination();
        }
        fakeRun = false;
        //end BeamOn(0);
        
        SetRunIDCounter(n_event);
        
    }
    
    
    
    // this function is a protected member of G4MTRunManager but we need to access it
    // from Mu2eG4_module, so we must make it public here
    void Mu2eG4MTRunManager::Mu2eG4InitializeEventLoop(G4int n_event)
    {
        //MTkernel->SetUpDecayChannels();
        G4MTRunManager::GetMTMasterRunManagerKernel()->SetUpDecayChannels();
        
        numberOfEventToBeProcessed = n_event;
        numberOfEventProcessed = 0;
        
        my_nworkers = G4MTRunManager::GetNumberOfThreads();
        
        if(!fakeRun)
        {
            nSeedsUsed = 0;
            nSeedsFilled = 0;
            
            //if(verboseLevel>0)
            //{ timer->Start(); }
            
            n_select_msg = -1;
            selectMacro = "";
            
            //initialize seeds
            //If user did not implement InitializeSeeds,
            // use default: nSeedsPerEvent seeds per event
            if( eventModuloDef > 0 )
                {
                    eventModulo = eventModuloDef;
                    if(eventModulo > numberOfEventToBeProcessed/my_nworkers)
                    {
                        eventModulo = numberOfEventToBeProcessed/my_nworkers;
                        if(eventModulo<1) eventModulo =1;
                        G4ExceptionDescription msgd;
                        msgd << "Event modulo is reduced to " << eventModulo
                             << " to distribute events to all threads.";
                                G4Exception("G4MTRunManager::InitializeEventLoop()",
                                            "Run10035", JustWarning, msgd);
                    }
                }
            else
            {
                eventModulo = int(std::sqrt(double(numberOfEventToBeProcessed/my_nworkers)));
                if(eventModulo<1) eventModulo =1;
            }
            if ( InitializeSeeds(n_event) == false && n_event>0 )
            {
                G4RNGHelper* helper = G4RNGHelper::GetInstance();
                switch(seedOncePerCommunication)
                {
                    case 0:
                        nSeedsFilled = n_event;
                        break;
                    case 1:
                        nSeedsFilled = my_nworkers;
                        break;
                    case 2:
                        nSeedsFilled = n_event/eventModulo + 1;
                        break;
                    default:
                        G4ExceptionDescription msgd;
                        msgd << "Parameter value <" << seedOncePerCommunication
                             << "> of seedOncePerCommunication is invalid. It is reset to 0." ;
                        G4Exception("G4MTRunManager::InitializeEventLoop()",
                                    "Run10036", JustWarning, msgd);
                        seedOncePerCommunication = 0;
                        nSeedsFilled = n_event;
                }

                // Generates up to nSeedsMax seed pairs only.
                if(nSeedsFilled>nSeedsMax) nSeedsFilled=nSeedsMax;
                
                CLHEP::HepRandomEngine* mrnge = const_cast<CLHEP::HepRandomEngine*> (G4MTRunManager::GetMasterRunManager()->getMasterRandomEngine());
                
                mrnge->flatArray(nSeedsPerEvent*nSeedsFilled,randDbl);
                helper->Fill(randDbl,nSeedsFilled,n_event,nSeedsPerEvent);
                
            }
        }

        //Now initialize workers. Check if user defined a WorkerThreadInitialization
        if ( userWorkerThreadInitialization == 0 )
        { userWorkerThreadInitialization = new G4UserWorkerThreadInitialization(); }
        
        //Prepare UI commands for threads
        PrepareCommandsStack();
        
        //Start worker threads
        CreateAndStartWorkers();
        
        // We need a barrier here. Wait for workers to start event loop.
        //This will return only when all workers have started processing events.
        WaitForReadyWorkers();
    }
    
    
    

    //this function is a protected member of G4MTRunManager but we need to access it
    //from Mu2eG4_module, so we must make it public here
    void Mu2eG4MTRunManager::Mu2eG4WaitForEndEventLoopWorkers()
    {
      if ( verboseLevel > 1 ) {
        G4cout << __func__ << " called" << G4endl;
      }
        WaitForEndEventLoopWorkers();
    }

    //we need control of the event loop in order to correctly break up the stages
    //to fit within the art framework.  we were getting a hang using the G4MTRunManager
    //version of RunTermination due to the call of WaitForEndEventLoopWorkers()
    void Mu2eG4MTRunManager::Mu2eG4RunTermination()
    {

      if ( verboseLevel > 0 ) {
        G4cout << __func__ << " called" << G4endl;
      }

        G4RunManager::TerminateEventLoop();
        G4RunManager::RunTermination();
    }
    
    void Mu2eG4MTRunManager::TestFunc()
    {
        
        //std::stringstream ss;
        //ss << std::this_thread::get_id();
        //int id = std::stoi(ss.str());
        
        std::cout << "Calling Mu2eG4MTRunManager::TestFunc()" << std::endl;
        
    }
    
    
} // end namespace mu2e
