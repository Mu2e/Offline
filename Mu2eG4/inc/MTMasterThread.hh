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


//G4 includes
//#include "G4UserWorkerInitialization.hh"

//C++ includes
#include <iostream>
#include <thread>
#include <pthread.h>
#include <mutex>
#include <condition_variable>

namespace fhicl { class ParameterSet; }

namespace mu2e {
    
    class Mu2eG4MTRunManager;

class MTMasterThread
{
  public:
    MTMasterThread(const fhicl::ParameterSet& pset);
    ~MTMasterThread();

    
    void beginRun() const;
    void endRun() const;
    void stopThread();
    
    inline Mu2eG4MTRunManager& masterRunManager() const { return *m_masterRunManager; }
    inline Mu2eG4MTRunManager* masterRunManagerPtr() const { return m_masterRunManager.get(); }
    
private:
    
    const fhicl::ParameterSet& pset_;
    
    void readES() const;
    
    enum class ThreadState { NotExist = 0, BeginRun = 1, EndRun = 2, Destruct = 3 };
    
    //const bool m_pUseMagneticField;
    //const bool m_pGeoFromDD4hep;
    
    std::shared_ptr<Mu2eG4MTRunManager> m_masterRunManager;
    std::thread m_masterThread;
    
    // ES products needed for Geant4 initialization
    /*mutable edm::ESWatcher<IdealGeometryRecord> idealGeomRcdWatcher_;
    mutable edm::ESWatcher<IdealMagneticFieldRecord> idealMagRcdWatcher_;
    mutable const DDCompactView* m_pDD;
    mutable const cms::DDCompactView* m_pDD4hep;
    mutable const MagneticField* m_pMF;
    mutable const HepPDT::ParticleDataTable* m_pTable;
    */
    mutable std::mutex m_protectMutex;
    mutable std::mutex m_threadMutex;
    mutable std::condition_variable m_notifyMasterCV;
    mutable std::condition_variable m_notifyMainCV;
    
    mutable ThreadState m_masterThreadState;
    
    mutable bool m_masterCanProceed;
    mutable bool m_mainCanProceed;
    mutable bool m_firstRun;
    mutable bool m_stopped;
    
};


}  // end namespace mu2e
#endif /* MT_MasterThread_hh */
