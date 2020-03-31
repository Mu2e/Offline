#ifndef MT_MasterThread_hh
#define MT_MasterThread_hh
//
// MTMasterThread.hh provides declarations for MTMasterThread
// class for the Mu2e G4 simulation.  MTMasterThread controls
// the workflow when running in MT mode.
//
// Author: Lisa Goodeough
// Date: 2019/07/17
//
//


//C++ includes

#include <thread>
#include <mutex>
#include <condition_variable>

#include "Mu2eG4/inc/Mu2eG4Config.hh"

namespace art { class Run; }

namespace mu2e {

  class Mu2eG4MTRunManager;
  class PhysicalVolumeHelper;

  class MTMasterThread
  {
  public:
    explicit MTMasterThread(const Mu2eG4Config::Top& conf);
    ~MTMasterThread();


    void beginRun() const;
    void endRun() const;
    void stopThread();
    void storeRunNumber(int art_runnumber);
    void readRunData(PhysicalVolumeHelper* phys_vol_help) const;

    inline Mu2eG4MTRunManager& masterRunManager() const { return *m_masterRunManager; }
    inline Mu2eG4MTRunManager* masterRunManagerPtr() const { return m_masterRunManager.get(); }

  private:

    bool m_mtDebugOutput;

    enum class ThreadState { NotExist = 0, BeginRun = 1, EndRun = 2, Destruct = 3 };

    std::shared_ptr<Mu2eG4MTRunManager> m_masterRunManager;
    std::thread m_masterThread;

    mutable std::mutex m_protectMutex;
    mutable std::mutex m_threadMutex;
    mutable std::condition_variable m_notifyMasterCV;
    mutable std::condition_variable m_notifyMainCV;

    mutable ThreadState m_masterThreadState;

    mutable bool m_masterCanProceed;
    mutable bool m_mainCanProceed;
    mutable bool m_firstRun;
    mutable bool m_stopped;

    int run_number;

  };


}  // end namespace mu2e
#endif /* MT_MasterThread_hh */
