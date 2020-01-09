//
//
// Author: Lisa Goodeough
// Date: 2019/07/17
//
//

//Mu2e includes
#include "Mu2eG4/inc/MTMasterThread.hh"
#include "Mu2eG4/inc/Mu2eG4MTRunManager.hh"

//art includes
#include "art/Framework/Principal/Run.h"
#include "fhiclcpp/ParameterSet.h"

//C++ includes
#include <iostream>

using namespace std;

namespace mu2e {

MTMasterThread::MTMasterThread(const fhicl::ParameterSet& pset)
    :
    pset_(pset),
    m_masterThreadState(ThreadState::NotExist),
    m_masterCanProceed(false),
    m_mainCanProceed(false),
    m_firstRun(true),
    m_stopped(false),
    run_number(0)
    {
        // Lock the mutex
        std::unique_lock<std::mutex> lk(m_threadMutex);
        
        // Create Genat4 master thread using a Lambda expression
        m_masterThread = std::thread([&]() {
           
            // Initialization
            std::shared_ptr<Mu2eG4MTRunManager> masterRunManager;
            
            // Lock the mutex (i.e. wait until the creating thread has called cv.wait()
            std::unique_lock<std::mutex> lk2(m_threadMutex);
            
            // Create the master run manager, and share it with the art::mu2e thread
            masterRunManager = std::make_shared<Mu2eG4MTRunManager>(pset);
            m_masterRunManager = masterRunManager;
            
            // State loop
            bool isG4Alive = false;
            while (true) {
                // Signal main thread that it can proceed
                m_mainCanProceed = true;
                std::cout << "Master thread: State loop, notify main thread" << std::endl;
                m_notifyMainCV.notify_one();
                
                // Wait until the main thread sends signal
                m_masterCanProceed = false;
                std::cout << "Master thread: State loop, starting wait" << std::endl;
                m_notifyMasterCV.wait(lk2, [&] { return m_masterCanProceed; });
                
                // Act according to the state
                std::cout << "Master thread: Woke up, state is " << static_cast<int>(m_masterThreadState) << std::endl;
                
                if (m_masterThreadState == ThreadState::BeginRun) {
                    // Initialize Geant4
                    std::cout << "Master thread: Initializing Geant4" << std::endl;
                    masterRunManager->initializeG4(run_number);
                    isG4Alive = true;
                } else if (m_masterThreadState == ThreadState::EndRun) {
                    // Stop Geant4
                    std::cout << "Master thread: Stopping Geant4" << std::endl;
                    masterRunManager->stopG4();
                    isG4Alive = false;
                } else if (m_masterThreadState == ThreadState::Destruct) {
                    std::cout << "Master thread: Breaking out of state loop" << std::endl;
                    if (isG4Alive)
                        throw cet::exception("G4STATELOGICERROR")
                        << "Geant4 is still alive, master thread state must be set to EndRun before Destruct\n";
                    break;
                } else {
                    throw cet::exception("G4STATELOGICERROR")
                    << "MTMasterThread: Illegal master thread state " << static_cast<int>(m_masterThreadState);
                }
            }
            
            // Cleanup
            std::cout << "MTMasterThread: start Mu2eG4MTRunManager destruction\n";
            std::cout << "Master thread: Am I unique owner of masterRunManager? "
                      << masterRunManager.unique() << std::endl;
            
            masterRunManager.reset();
            //G4PhysicalVolumeStore::Clean();
            
            std::cout << "Master thread: reset shared_ptr" << std::endl;
            lk2.unlock();
            std::cout << "MTMasterThread: Master thread is finished" << std::endl;
        });
        
        // Start waiting for a signal from the condition variable (releases the lock temporarily)
        // First for initialization
        m_mainCanProceed = false;
        std::cout << "Main thread: Signal master for initialization" << std::endl;
        m_notifyMainCV.wait(lk, [&]() { return m_mainCanProceed; });
        
        lk.unlock();
        std::cout << "MTMasterThread: Master thread is constructed" << std::endl;
    }
    

MTMasterThread::~MTMasterThread()
    {
        if (!m_stopped) {
            stopThread();
        }
    }
 
    
void MTMasterThread::beginRun() const {
    
        std::lock_guard<std::mutex> lk(m_protectMutex);
        std::unique_lock<std::mutex> lk2(m_threadMutex);
    
        m_masterThreadState = ThreadState::BeginRun;
        m_masterCanProceed = true;
        m_mainCanProceed = false;
        std::cout << "MTMasterThread: Signal master for BeginRun" << std::endl;
        m_notifyMasterCV.notify_one();
        m_notifyMainCV.wait(lk2, [&]() { return m_mainCanProceed; });
        
        lk2.unlock();
        std::cout << "MTMasterThread: finish BeginRun" <<std::endl;
    }
 

void MTMasterThread::endRun() const {
        
        std::lock_guard<std::mutex> lk(m_protectMutex);
        std::unique_lock<std::mutex> lk2(m_threadMutex);
        
        m_masterThreadState = ThreadState::EndRun;
        m_mainCanProceed = false;
        m_masterCanProceed = true;
        std::cout << "MTMasterThread: Signal master for EndRun" <<std::endl;
        m_notifyMasterCV.notify_one();
        m_notifyMainCV.wait(lk2, [&]() { return m_mainCanProceed; });
        lk2.unlock();
        std::cout << "MTMasterThread: finish EndRun" << std::endl;
    }

    
void MTMasterThread::stopThread() {
        if (m_stopped) {
            return;
        }
        std::cout << "MTMasterThread::stopThread: stop main thread" << std::endl;
        
        // Release our instance of the shared master run manager, so that
        // the G4 master thread can do the cleanup. Then notify the master
        // thread, and join it.
        std::unique_lock<std::mutex> lk2(m_threadMutex);
        m_masterRunManager.reset();
        std::cout << "Main thread: reset shared_ptr" << std::endl;
        
        m_masterThreadState = ThreadState::Destruct;
        m_masterCanProceed = true;
        std::cout << "MTMasterThread::stopThread: notify" << std::endl;
        m_notifyMasterCV.notify_one();
        lk2.unlock();
        
        std::cout << "Main thread: joining master thread" << std::endl;
        m_masterThread.join();
        std::cout << "MTMasterThread::stopThread: main thread finished" << std::endl;
        m_stopped = true;
    }
    
  
void MTMasterThread::storeRunNumber(int art_runnumber)
    {
        run_number = art_runnumber;
    }
 
    
void MTMasterThread::readRunData(PhysicalVolumeHelper* phys_vol_help) const {
    
    m_masterRunManager->setPhysVolumeHelper(phys_vol_help);
    
}
    
}// end namespace mu2e



