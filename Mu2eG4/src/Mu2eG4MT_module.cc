// A Producer Module that runs Geant4 and adds its output to the event.
//
// Original author Rob Kutschke
//
// Notes:



// Mu2e includes
#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Offline/Mu2eG4/inc/exportG4PDT.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"
#include "Offline/Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4Inputs.hh"
#include "Offline/Mu2eG4/inc/checkConfigRelics.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4WorkerRunManager.hh"
#include "Offline/Mu2eG4/inc/MTMasterThread.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4IOConfigHelper.hh"
#include "Offline/Mu2eG4/inc/writePhysicalVolumes.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4MTRunManager.hh"
#include "Offline/Mu2eG4/inc/validGeometryOrThrow.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ScoringManager.hh"
#include "Offline/SeedService/inc/SeedService.hh"

// Data products that will be produced by this module.
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Offline/MCDataProducts/inc/StatusG4.hh"

// From art and its tool chain.
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Utilities/Globals.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Geant4 includes
#include "Geant4/G4Run.hh"
#include "Geant4/G4VUserPhysicsList.hh"
#include "Geant4/G4ScoringManager.hh"

// C++ includes.
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

// TBB includes
#include "tbb/concurrent_hash_map.h"
#include "tbb/task_group.h"

using namespace std;

namespace {
  art::InputTag const invalid_tag{};
}

namespace tbb {
  template<>
  struct tbb::detail::d1::tbb_hash_compare<std::thread::id> {
    tbb_hash_compare() {}
    static size_t hash(std::thread::id tid) {
      std::ostringstream oss;
      oss << tid;
      return std::stoull(oss.str());
    }
    static bool equal(std::thread::id const k1, std::thread::id const k2) {
      return k1==k2;
    }
  };
}

namespace mu2e {

  class Mu2eG4MT : public art::SharedProducer {
  public:
    using Parameters = art::SharedProducer::Table<Mu2eG4Config::Top>;

    Mu2eG4MT(Parameters const& pars, art::ProcessingFrame const& pf);
    // Accept compiler supplied d'tor

  private:
    void produce(art::Event& e, art::ProcessingFrame const& pf) override;
    void endJob(art::ProcessingFrame const& pf) override;
    void beginRun(art::Run &r, art::ProcessingFrame const& pf) override;
    void endRun(art::Run &r, art::ProcessingFrame const& pf) override;
    void beginSubRun(art::SubRun &sr, art::ProcessingFrame const& pf) override;
    void endSubRun(art::SubRun &sr, art::ProcessingFrame const& pf) override;

    Mu2eG4Config::Top conf_;
    Mu2eG4ResourceLimits mu2elimits_;
    Mu2eG4TrajectoryControl trajectoryControl_;
    Mu2eG4Inputs multiStagePars_;

    unsigned simStage_;

    std::unique_ptr<MTMasterThread> masterThread;
    std::unique_ptr<Mu2eG4ScoringManager> _scorer;

    // Do we issue warnings about multiple runs?
    bool _warnEveryNewRun;

    // Do we want to export the G4 particle data table.
    bool  _exportPDTStart;
    bool  _exportPDTEnd;

    std::string storePhysicsTablesDir_;

    int _rmvlevel;
    int _mtDebugOutput;

    // Instance name of the timeVD StepPointMC data product.
    const StepInstanceName _tvdOutputName;
    bool timeVD_enabled_;

    // Helps with indexology related to persisting G4 volume information.
    // string to ptr maps, speed optimization
    // if in MT mode, only allow lookup, don't allow add
    // do a counter that counts how mnay times it was called with an unknown process
    //NOW DONE IN Master Run Manager code
    PhysicalVolumeHelper physVolHelper_;
    Mu2eG4IOConfigHelper ioconf_;

    // Do the G4 initialization that must be done only once per job, not once per run
    void initializeG4( GeometryService& geom, art::Run const& run );

    const bool standardMu2eDetector_;
    G4ThreeVector originInWorld;

    // count the number of events that have been excluded because they did not
    // pass the filtering in Mu2eG4EventAction
    int numExcludedEvents = 0;

    int const num_schedules{art::Globals::instance()->nschedules()};
    int const num_threads{art::Globals::instance()->nthreads()};

    typedef tbb::concurrent_hash_map< std::thread::id, std::unique_ptr<Mu2eG4WorkerRunManager> > WorkerRMMap;
    WorkerRMMap myworkerRunManagerMap;
  }; // end G4 header


  Mu2eG4MT::Mu2eG4MT(Parameters const& pars, art::ProcessingFrame const& procFrame):
    SharedProducer{pars},
    conf_(pars()),
    mu2elimits_(pars().ResourceLimits()),
    trajectoryControl_(pars().TrajectoryControl()),
    multiStagePars_(pars().inputs()),
    simStage_(-1u),

    masterThread(std::make_unique<MTMasterThread>(pars(),mu2elimits_ )),
    _scorer(std::make_unique<Mu2eG4ScoringManager>(G4ScoringManager::GetScoringManager(),
                                                   conf_.scoring(),conf_.physics(),conf_.debug())),

    _warnEveryNewRun(pars().debug().warnEveryNewRun()),
    _exportPDTStart(pars().debug().exportPDTStart()),
    _exportPDTEnd(pars().debug().exportPDTEnd()),

    storePhysicsTablesDir_(pars().debug().storePhysicsTablesDir()),

    _rmvlevel(pars().debug().diagLevel()),
    _mtDebugOutput(pars().debug().mtDebugOutput()),

    _tvdOutputName(StepInstanceName::timeVD),
    timeVD_enabled_(pars().SDConfig().TimeVD().enabled()),
    physVolHelper_(),
    ioconf_(pars(), producesCollector(), consumesCollector()),
    standardMu2eDetector_((art::ServiceHandle<GeometryService>())->isStandardMu2eDetector())
    {
      // produces() and consumes()  calls are handled by Mu2eG4IOConfigHelper

      // This statement requires that the external libraries the module uses are thread-safe,
      // and that the data member members are used in a thread-safe manner
      async<art::InEvent>();

      if (num_schedules>1) {
        G4cout << "Mu2eG4MT starting "<< num_schedules <<" threads" <<endl;
      }
    } // end Mu2eG4MT constructor



  void Mu2eG4MT::beginRun( art::Run &run, art::ProcessingFrame const& procFrame) {

    art::ServiceHandle<GeometryService> geom;
    SimpleConfig const& config  = geom->config();
    checkConfigRelics(config);

    static int ncalls(0);
    ++ncalls;

    // Do the main initialization of G4; only once per job.
    if ( ncalls == 1 ) {
      initializeG4( *geom, run );
    } else {
      if ( ncalls ==2 || _warnEveryNewRun ){
        mf::LogWarning log("G4");
        log << "G4 does not change state when we cross run boundaries - hope this is OK .... ";
        if ( ncalls == 2 && !_warnEveryNewRun ){
          log << "\nThis message will not be repeated on subsequent new runs.";
        }
      }
    }

    // Some end-of-begin run things that need to be done only once per job.
    if ( ncalls == 1 ){
      validGeometryOrThrow( _rmvlevel );
      if ( _exportPDTStart) {
        exportG4PDT( "Start:" );
      }
    }

  }//Mu2eG4MT::beginRun



  void Mu2eG4MT::initializeG4( GeometryService& geom, art::Run const& run ) {
    if (standardMu2eDetector_) {
      geom.addWorldG4(*GeomHandle<Mu2eHall>());
      originInWorld = GeomHandle<WorldG4>()->mu2eOriginInWorld();
    }

    if ( _rmvlevel > 0 ) {
      mf::LogInfo logInfo("GEOM");
      logInfo << "Initializing Geant4 for " << run.id()
              << " with verbosity " << _rmvlevel << endl;
    }

    _scorer->initialize();

    masterThread->storeRunNumber(run.id().run());
    masterThread->readRunData(&physVolHelper_);
    masterThread->beginRun();

  }//Mu2eG4MT::initializeG4


  void Mu2eG4MT::beginSubRun(art::SubRun& sr, art::ProcessingFrame const& procFrame) {
    if(multiStagePars_.simStageOverride()) {
      simStage_ = *multiStagePars_.simStageOverride();
    }
    else {
      std::optional<art::InputTag> in;
      if(multiStagePars_.multiStage()) {
        in.emplace(multiStagePars_.inputPhysVolumeMultiInfo());
      }
      simStage_ = writePhysicalVolumes(sr,
                                       in,
                                       physVolHelper_.persistentSingleStageInfo(),
                                       "");
    }
  }

  void Mu2eG4MT::endSubRun(art::SubRun& sr, art::ProcessingFrame const& procFrame) {
    if(multiStagePars_.simStageOverride()) {
      const unsigned pvstage =
        writePhysicalVolumes(sr,
                             multiStagePars_.inputPhysVolumeMultiInfo(),
                             physVolHelper_.persistentSingleStageInfo(),
                             "");

      if(pvstage != simStage_) {
        throw cet::exception("BADINPUT")
          << "Mu2eG4MT::endSubRun() Error: inconsistent simStage: "
          <<simStage_<<" vs "<<pvstage<<"\n";
      }
    }
   _scorer->dumpInDataProduct(sr);
   _scorer->reset();
  }

  // Create one G4 event and copy its output to the art::event.
  void Mu2eG4MT::produce(art::Event& event, art::ProcessingFrame const& procFrame) {

    int schedID = std::stoi(std::to_string(procFrame.scheduleID().id()));
    auto const tid = std::this_thread::get_id();

    WorkerRMMap::accessor access_workerMap;

    if (!myworkerRunManagerMap.find(access_workerMap, tid)){
      if (_mtDebugOutput > 0){
        G4cout << "FOR TID: " << tid << ", NO WORKER.  We are making one.\n";
      }
      myworkerRunManagerMap.insert(access_workerMap, tid);
      access_workerMap->second = std::make_unique<Mu2eG4WorkerRunManager>(conf_, ioconf_, tid);
    }

    if (event.id().event() == 1 && _mtDebugOutput > 0) {
      G4cout << "Our RMmap has " << myworkerRunManagerMap.size() << " members\n";
    }

    myworkerRunManagerMap.find(access_workerMap, tid);
    Mu2eG4WorkerRunManager* scheduleWorkerRM = (access_workerMap->second).get();
    access_workerMap.release();

    if (_mtDebugOutput > 1){
      G4cout << "FOR SchedID: " << schedID << ", TID=" << tid << ", workerRunManagers[schedID].get() is:" << scheduleWorkerRM << "\n";
    }

    //if this is the first time the thread is being used, it should be initialized
    if (!scheduleWorkerRM->workerRMInitialized()){
      scheduleWorkerRM->initializeThread(masterThread->masterRunManagerPtr(), originInWorld);
      scheduleWorkerRM->initializeRun(&event);
    }

    Mu2eG4PerThreadStorage* perThreadStore = scheduleWorkerRM->getMu2eG4PerThreadStorage();
    perThreadStore->initializeEventInfo(&event, simStage_);
    scheduleWorkerRM->processEvent(event.id());

    if (_mtDebugOutput > 2){
      G4cout << "Current Event in RM is: " << scheduleWorkerRM->GetCurrentEvent()->GetEventID() << "\n";
    }

    /////////////////////////////////////////////////////////////////////////////////////
    //putting data into the event

    if(!perThreadStore->eventPassed()) {
      perThreadStore->clearData();
      numExcludedEvents++;
    }
    else {
      perThreadStore->putDataIntoEvent();

      if(multiStagePars_.updateEventLevelVolumeInfos()) {
        const unsigned pvstage =
          writePhysicalVolumes(event,
                               multiStagePars_.updateEventLevelVolumeInfos()->input,
                               physVolHelper_.persistentSingleStageInfo(),
                               multiStagePars_.updateEventLevelVolumeInfos()->outInstance);

        if(pvstage != simStage_) {
          throw cet::exception("BADINPUT")
            << "Mu2eG4MT::produce() Error: inconsistent simStage: "
            <<simStage_<<" vs "<<pvstage<<"\n";
        }
      }
    }

    scheduleWorkerRM->TerminateOneEvent();

  }//end Mu2eG4MT::produce


  // Tell G4 that this run is over.
  void Mu2eG4MT::endRun(art::Run & run, art::ProcessingFrame const& procFrame) {

    if (_mtDebugOutput > 0){
      G4cout << "At endRun pt1, we have " << myworkerRunManagerMap.size() << " members in the map "
             << "and are running " << num_threads << " threads.\n" ;
    }
    else if (num_threads < static_cast <int> (myworkerRunManagerMap.size()) && _mtDebugOutput > 0){
      G4cout << "At endRun pt1, we have " << myworkerRunManagerMap.size() << " members in the map "
             << "and are running " << num_threads << " threads.\n" ;
    }


    if (storePhysicsTablesDir_!="") {
      if ( _rmvlevel > 0 ) {
        G4cout << __func__ << " Will write out physics tables to "
               << storePhysicsTablesDir_
               << G4endl;
      }
      masterThread->masterRunManagerPtr()->getMasterPhysicsList()->StorePhysicsTable(storePhysicsTablesDir_);
    }

    std::atomic<int> threads_left = num_threads;
    tbb::task_group g;
    for (int i = 0; i < num_threads; ++i) {

      auto destroy_worker = [&threads_left, i, this] {
        WorkerRMMap::accessor access_workerMap;
        std::thread::id this_tid = std::this_thread::get_id();

        if (myworkerRunManagerMap.find(access_workerMap, this_tid)) {
          access_workerMap->second.reset();
          myworkerRunManagerMap.erase(access_workerMap);
        }

        access_workerMap.release();
        --threads_left;
        while (threads_left != 0) {}
        return;
      };
      g.run(destroy_worker);
    }//for
    g.wait();

    if (_mtDebugOutput > 0){
      G4cout << "At endRun pt2, we have " << myworkerRunManagerMap.size() << " members in the map.\n";
    }

    //This cleans up the worker run managers that are in threads no longer being used, i.e. 'transient threads'
    WorkerRMMap::iterator it = myworkerRunManagerMap.begin();
    while (it != myworkerRunManagerMap.end()) {
      if (_mtDebugOutput > 0){
        G4cout << "releasing RM for thread ID" << it->first << std::endl;
      }
      it->second.release();
      ++it;
    }

    if ( _rmvlevel > 0 ) {
      G4cout << "at endRun: numExcludedEvents = " << numExcludedEvents << G4endl;
    }
    myworkerRunManagerMap.clear();
    masterThread->endRun();
  }



  void Mu2eG4MT::endJob(art::ProcessingFrame const& procFrame) {

    if ( _exportPDTEnd ) exportG4PDT( "End:" );
    physVolHelper_.endRun();
  }


} // End of namespace mu2e

DEFINE_ART_MODULE(mu2e::Mu2eG4MT)
