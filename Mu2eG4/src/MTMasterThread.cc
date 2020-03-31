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

//C++ includes
#include <iostream>

using namespace std;

namespace mu2e {

  MTMasterThread::MTMasterThread(const Mu2eG4Config::Top& conf)
    :
    m_mtDebugOutput(conf.debug().mtDebugOutput()),
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
        masterRunManager = std::make_shared<Mu2eG4MTRunManager>(conf);
        m_masterRunManager = masterRunManager;

        // State loop
        bool isG4Alive = false;
        while (true) {
          // Signal main thread that it can proceed
          m_mainCanProceed = true;
          if (m_mtDebugOutput) {
            G4cout << "Master thread: State loop, notify main thread" << G4endl;
          }

          m_notifyMainCV.notify_one();

          // Wait until the main thread sends signal
          m_masterCanProceed = false;
          if (m_mtDebugOutput) {
            G4cout << "Master thread: State loop, starting wait" << G4endl;
          }
          m_notifyMasterCV.wait(lk2, [&] { return m_masterCanProceed; });

          // Act according to the state
          if (m_mtDebugOutput) {
            G4cout << "Master thread: Woke up, state is " << static_cast<int>(m_masterThreadState) << G4endl;
          }

          if (m_masterThreadState == ThreadState::BeginRun) {
            // Initialize Geant4
            if (m_mtDebugOutput) {
              G4cout << "Master thread: Initializing Geant4" << G4endl;
            }
            masterRunManager->initializeG4(run_number);
            isG4Alive = true;
          } else if (m_masterThreadState == ThreadState::EndRun) {
            // Stop Geant4
            if (m_mtDebugOutput) {
              G4cout << "Master thread: Stopping Geant4" << G4endl;
            }
            masterRunManager->stopG4();
            isG4Alive = false;
          } else if (m_masterThreadState == ThreadState::Destruct) {
            if (m_mtDebugOutput) {
              G4cout << "Master thread: Breaking out of state loop" << G4endl;
            }
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
        if (m_mtDebugOutput) {
          G4cout << "MTMasterThread: start Mu2eG4MTRunManager destruction\n";
          G4cout << "Master thread: Am I unique owner of masterRunManager? "
          << masterRunManager.unique() << G4endl;
        }

        masterRunManager.reset();
        //G4PhysicalVolumeStore::Clean();

        if (m_mtDebugOutput) {
          G4cout << "Master thread: reset shared_ptr" << G4endl;
        }
        lk2.unlock();
        if (m_mtDebugOutput) {
          G4cout << "MTMasterThread: Master thread is finished" << G4endl;
        }
      });

    // Start waiting for a signal from the condition variable (releases the lock temporarily)
    // First for initialization
    m_mainCanProceed = false;
    if (m_mtDebugOutput) {
      G4cout << "Main thread: Signal master for initialization" << G4endl;
    }
    m_notifyMainCV.wait(lk, [&]() { return m_mainCanProceed; });

    lk.unlock();
    if (m_mtDebugOutput) {
      G4cout << "MTMasterThread: Master thread is constructed" << G4endl;
    }
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
    if (m_mtDebugOutput) {
      G4cout << "MTMasterThread: Signal master for BeginRun" << G4endl;
    }
    m_notifyMasterCV.notify_one();
    m_notifyMainCV.wait(lk2, [&]() { return m_mainCanProceed; });

    lk2.unlock();
    if (m_mtDebugOutput) {
      G4cout << "MTMasterThread: finish BeginRun" << G4endl;
    }
  }


  void MTMasterThread::endRun() const {

    std::lock_guard<std::mutex> lk(m_protectMutex);
    std::unique_lock<std::mutex> lk2(m_threadMutex);

    m_masterThreadState = ThreadState::EndRun;
    m_mainCanProceed = false;
    m_masterCanProceed = true;
    if (m_mtDebugOutput) {
      G4cout << "MTMasterThread: Signal master for EndRun" <<G4endl;
    }
    m_notifyMasterCV.notify_one();
    m_notifyMainCV.wait(lk2, [&]() { return m_mainCanProceed; });
    lk2.unlock();
    if (m_mtDebugOutput) {
      G4cout << "MTMasterThread: finish EndRun" << G4endl;
    }
  }


  void MTMasterThread::stopThread() {
    if (m_stopped) {
      return;
    }
    if (m_mtDebugOutput) {
      G4cout << "MTMasterThread::stopThread: stop main thread" << G4endl;
    }

    // Release our instance of the shared master run manager, so that
    // the G4 master thread can do the cleanup. Then notify the master
    // thread, and join it.
    std::unique_lock<std::mutex> lk2(m_threadMutex);
    m_masterRunManager.reset();
    if (m_mtDebugOutput) {
      G4cout << "Main thread: reset shared_ptr" << G4endl;
    }

    m_masterThreadState = ThreadState::Destruct;
    m_masterCanProceed = true;
    if (m_mtDebugOutput) {
      G4cout << "MTMasterThread::stopThread: notify" << G4endl;
    }
    m_notifyMasterCV.notify_one();
    lk2.unlock();

    if (m_mtDebugOutput) {
      G4cout << "Main thread: joining master thread" << G4endl;
    }
    m_masterThread.join();
    if (m_mtDebugOutput) {
      G4cout << "MTMasterThread::stopThread: main thread finished" << G4endl;
    }
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
