// A Producer Module that runs Geant4 and adds its output to the event.
//
// Original author Rob Kutschke
//
// Notes:



// Mu2e includes
#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/Mu2eStudyWorld.hh"
#include "Mu2eG4/inc/IMu2eG4Cut.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/exportG4PDT.hh"

#include "Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

//#include "GeometryService/inc/WorldG4.hh"
//#include "Mu2eG4/inc/ActionInitialization.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/physicsListDecider.hh"
//#include "Mu2eG4/inc/preG4InitializeTasks.hh"

#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

//#include "Mu2eG4/inc/ExtMonFNALPixelSD.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eG4/inc/generateFieldMap.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"
#include "Mu2eG4/inc/Mu2eG4MultiStageParameters.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/checkConfigRelics.hh"
#include "Mu2eG4/inc/Mu2eG4MTRunManager.hh"
#include "Mu2eG4/inc/Mu2eG4WorkerRunManager.hh"
#include "Mu2eG4/inc/MTMasterThread.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"



// Data products that will be produced by this module.
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"

// From art and its tool chain.
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Utilities/Globals.h"

// Geant4 includes

#include "G4Run.hh"
#include "G4Timer.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ParticleHPManager.hh"
#include "G4HadronicProcessStore.hh"
#include "G4SDManager.hh"


// C++ includes.
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <utility>
#include <mutex>

// TBB includes
#include "tbb/concurrent_hash_map.h"

using namespace std;

namespace {
  art::InputTag const invalid_tag{};
}

namespace mu2e {

  class Mu2eG4 : public art::SharedProducer {
  public:
    Mu2eG4(fhicl::ParameterSet const& pSet, art::ProcessingFrame const& pf);
    // Accept compiler supplied d'tor

  private:
    void produce(art::Event& e, art::ProcessingFrame const& pf) override;
    void endJob(art::ProcessingFrame const& pf) override;
    void beginRun(art::Run &r, art::ProcessingFrame const& pf) override;
    void endRun(art::Run &r, art::ProcessingFrame const& pf) override;
    void beginSubRun(art::SubRun &sr, art::ProcessingFrame const& pf) override;

    fhicl::ParameterSet pset_;
    Mu2eG4ResourceLimits mu2elimits_;
    Mu2eG4MultiStageParameters multiStagePars_;

    // The THREE functions that call new G4RunManger functions and break G4's BeamOn() into 3 pieces
    //void BeamOnBeginRun( unsigned int runNumber);
    //void BeamOnDoOneArtEvent( int eventNumber, G4int, const char* macroFile=0, G4int n_select=-1 );
    void BeamOnEndRun();

      
    std::unique_ptr<MTMasterThread> masterThread;
    std::vector< std::unique_ptr<Mu2eG4WorkerRunManager> > workerRunManagers;
      
    // Do we issue warnings about multiple runs?
    bool _warnEveryNewRun;

    // Do we want to export the G4 particle data table.
    bool  _exportPDTStart;
    bool  _exportPDTEnd;
      

    //these cut objects are used in the master thread to indicate what data product is produced
    //additional thread-local cut objects are owned by ActionInitialization
    std::unique_ptr<IMu2eG4Cut> stackingCuts_;
    std::unique_ptr<IMu2eG4Cut> steppingCuts_;
    std::unique_ptr<IMu2eG4Cut> commonCuts_;

    int _rmvlevel;
    int _tmvlevel;
    int _smvlevel;

    art::InputTag _generatorModuleLabel;

    // Helps with indexology related to persisting G4 volume information.
    // string to ptr maps, speed optimization
    // if in MT mode, only allow lookup, don't allow add
    // do a counter that counts how mnay times it was called with an unknown process
    //NOW DONE IN Master Run Manager code
    PhysicalVolumeHelper physVolHelper_;

//    vector< SensitiveDetectorHelper > SensitiveDetectorHelpers;
//    ExtMonFNALPixelSD       *_extMonFNALPixelSD;

    SensitiveDetectorHelper sensitiveDetectorHelper_;
      
      
    // Do the G4 initialization that must be done only once per job, not once per run
    void initializeG4( GeometryService& geom, art::Run const& run );

    unique_ptr<G4Timer> _timer; // local Mu2e per Geant4 event timer
    // Counters for cumulative time spent processing events by Geant4
    G4double _realElapsed;
    G4double _systemElapsed;
    G4double _userElapsed;

    const bool standardMu2eDetector_;
    G4ThreeVector originInWorld;
      
    // count the number of events that have been excluded because they did not
    // pass the filtering in Mu2eG4EventAction
    int numExcludedEvents = 0;
    
    std::mutex initmutex_;
      
    CLHEP::HepJamesRandom _engine;

    int num_schedules;
    
    typedef tbb::concurrent_hash_map<int, Mu2eG4PerThreadStorage> ScheduleMap;
    ScheduleMap myPerSchedMap;
    ScheduleMap::accessor access_PerSchedMap;

  }; // end G4 header

    
  Mu2eG4::Mu2eG4(fhicl::ParameterSet const& pSet, art::ProcessingFrame const& procFrame):
    SharedProducer{pSet},
    pset_(pSet),
    mu2elimits_(pSet.get<fhicl::ParameterSet>("ResourceLimits")),
    multiStagePars_(pSet.get<fhicl::ParameterSet>("MultiStageParameters")),

    masterThread(std::make_unique<MTMasterThread>(pSet)),
    

    _warnEveryNewRun(pSet.get<bool>("debug.warnEveryNewRun",false)),
    _exportPDTStart(pSet.get<bool>("debug.exportPDTStart",false)),
    _exportPDTEnd(pSet.get<bool>("debug.exportPDTEnd",false)),

    stackingCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4StackingOnlyCut", {}), mu2elimits_)),
    steppingCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4SteppingOnlyCut", {}), mu2elimits_)),
    commonCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4CommonCut", {}), mu2elimits_)),

    _rmvlevel(pSet.get<int>("debug.diagLevel",0)),
    _tmvlevel(pSet.get<int>("debug.trackingVerbosityLevel",0)),
    _smvlevel(pSet.get<int>("debug.steppingVerbosityLevel",0)),
    
    _generatorModuleLabel(pSet.get<string>("generatorModuleLabel", "")),
    physVolHelper_(),
    sensitiveDetectorHelper_(pSet.get<fhicl::ParameterSet>("SDConfig", fhicl::ParameterSet())),
//    _extMonFNALPixelSD(),

    _timer(make_unique<G4Timer>()),
    _realElapsed(0.),
    _systemElapsed(0.),
    _userElapsed(0.),
    standardMu2eDetector_((art::ServiceHandle<GeometryService>())->isStandardMu2eDetector()),
    
    //NEED TO FIGURE OUT HOW TO CONNECT THIS ENGINE TO THE G4 ENGINE
    _engine{art::ServiceHandle<SeedService>{}->getSeed()}
    {
        if((_generatorModuleLabel == art::InputTag()) && multiStagePars_.genInputHits().empty()) {
            throw cet::exception("CONFIG")
            << "Error: both generatorModuleLabel and genInputHits are empty - nothing to do!\n";
        }
        
   //         _runManager.reset(new Mu2eG4MTRunManager());
 
    // This statement requires that the external libraries the module uses are thread-safe,
    // and that the data member members are used in a thread-safe manner
    async<art::InEvent>();
        
    auto& collector = producesCollector();
        
    sensitiveDetectorHelper_.declareProducts(collector);
        
    stackingCuts_->declareProducts(collector);
    steppingCuts_->declareProducts(collector);
    commonCuts_->declareProducts(collector);
       
    produces<StatusG4>();
    produces<SimParticleCollection>();
     
    if (_generatorModuleLabel != invalid_tag) {
      consumes<GenParticleCollection>(_generatorModuleLabel);
    }
    

    // The string "G4Engine" is magic; see the docs for RandomNumberGenerator.
    // This does not work in a Shared Module
    //createEngine( art::ServiceHandle<SeedService>()->getSeed(), "G4Engine");
        
        
    num_schedules = art::Globals::instance()->nschedules();
    workerRunManagers.resize(num_schedules);
    std::cout << "WE WILL RUN " << num_schedules << " SCHEDULES" <<  std::endl;
        
    for(int i=0; i<num_schedules; ++i)
    {
        myPerSchedMap.emplace(std::make_pair( int(i), Mu2eG4PerThreadStorage(pSet,i+10) ));
    }

} // end G4:G4(fhicl::ParameterSet const& pSet);


    
    
    
void Mu2eG4::beginRun( art::Run &run, art::ProcessingFrame const& procFrame){
    
   std::cout << "CALLING beginRun from SCHEDULE: " << procFrame.scheduleID() << std::endl;
    
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
    
    // A few more things that only need to be done only once per job,
    // not once per run, but which need to be done after the call to
    // BeamOnBeginRun.
    if ( ncalls == 1 ) {
        // Since the cuts that are put into the event are owned by the individual threads,
        // I do not think we need to do this.  Taking it out for now.
        //stackingCuts_->finishConstruction(originInWorld);
        //steppingCuts_->finishConstruction(originInWorld);
        //commonCuts_->finishConstruction(originInWorld);
        if ( _exportPDTStart ) exportG4PDT( "Start:" );//once per job
        }//if ( ncalls == 1 )
    
}//Mu2eG4::beginRun


    
    
    
    
void Mu2eG4::initializeG4( GeometryService& geom, art::Run const& run )
    {
        if (standardMu2eDetector_) {
            geom.addWorldG4(*GeomHandle<Mu2eHall>());
        }

        if ( _rmvlevel > 0 ) {
            mf::LogInfo logInfo("GEOM");
            logInfo << "Initializing Geant4 for " << run.id()
            << " with verbosity " << _rmvlevel << endl;
            logInfo << " Configured simParticleNumberOffset = "<< multiStagePars_.simParticleNumberOffset() << endl;
        }
        
        masterThread->storeRunNumber(run.id().run());
        masterThread->readRunData(&physVolHelper_);
        masterThread->beginRun();
        
    }//Mu2eG4::initializeG4

    
 


void Mu2eG4::beginSubRun(art::SubRun& sr, art::ProcessingFrame const& procFrame)
    {
/*        using Collection_t = PhysicalVolumeInfoMultiCollection;
        auto mvi = std::make_unique<Collection_t>();
        
        if(multiStagePars_.inputPhysVolumeMultiInfo() != invalid_tag) {
            // Copy over data from the previous simulation stages
            auto const& ih = sr.getValidHandle<Collection_t>(multiStagePars_.inputPhysVolumeMultiInfo());
            mvi->reserve(1 + ih->size());
            mvi->insert(mvi->begin(), ih->cbegin(), ih->cend());
            
        }
        
        // Append info for the current stage
        mvi->emplace_back(multiStagePars_.simParticleNumberOffset(), physVolHelper_.persistentSingleStageInfo());
        
        sr.put(std::move(mvi));
*/
}
        
        
        
// Create one G4 event and copy its output to the art::event.
void Mu2eG4::produce(art::Event& event, art::ProcessingFrame const& procFrame) {
    
    if (num_schedules>1) {
        if (   multiStagePars_.inputSimParticles() != art::InputTag()
                || multiStagePars_.inputMCTrajectories() != art::InputTag()
                || !(multiStagePars_.genInputHits().empty()) ) {
                throw cet::exception("CONFIG")
                << "Error: You are trying to run in MT mode with input from previous stages.  This is an invalid configuration!\n";
        }
    }
    
    art::Handle<GenParticleCollection> gensHandle;
    if(!(_generatorModuleLabel == art::InputTag())) {
        event.getByLabel(_generatorModuleLabel, gensHandle);
    }
    
    HitHandles genInputHits;
    for(const auto& i : multiStagePars_.genInputHits()) {
        genInputHits.emplace_back(event.getValidHandle<StepPointMCCollection>(i));
    }
    
    art::ProductID simPartId(event.getProductID<SimParticleCollection>());
    art::EDProductGetter const* simProductGetter = event.productGetter(simPartId);
    
    SimParticleHelper spHelper(multiStagePars_.simParticleNumberOffset(), simPartId, &event, simProductGetter);
    SimParticlePrimaryHelper parentHelper(&event, simPartId, gensHandle, simProductGetter);

    
    int schedID = std::stoi(std::to_string(procFrame.scheduleID().id()));
    
    if(workerRunManagers[schedID].get()==nullptr){
            std::cerr << "FOR SchedID: " << schedID << ", NO WORKER.  We are making one.\n";
            workerRunManagers[schedID].reset(new Mu2eG4WorkerRunManager(pset_));
    }

    Mu2eG4WorkerRunManager* scheduleWorkerRM = workerRunManagers[schedID].get();
    
    myPerSchedMap.find(access_PerSchedMap, schedID);
    Mu2eG4PerThreadStorage* perThreadStore = &(access_PerSchedMap->second);
    access_PerSchedMap.release();
    
        
    //if this is the first time the thread is being used, it should be initialized
    if (!scheduleWorkerRM->workerRMInitialized()){
        scheduleWorkerRM->initializePTS(perThreadStore);
        scheduleWorkerRM->initializeThread(masterThread->masterRunManagerPtr(), originInWorld);
        scheduleWorkerRM->initializeRun(&event);
    }
    
    perThreadStore->initializeEventInfo(&event, &spHelper, &parentHelper, &genInputHits, _generatorModuleLabel);
    scheduleWorkerRM->processEvent(&event);
        
        
/////////////////////////////////////////////////////////////////////////////////////
    std::cout << "Putting data into the EVENT" << std::endl;
    
    //std::unique_ptr<SimParticleCollection> simsToCheck = std::move(perThreadStore->getSimPartCollection());
    std::unique_ptr<SimParticleCollection> simsToCheck = perThreadStore->getSimPartCollection();
    
    if (simsToCheck == nullptr) {
        numExcludedEvents++;
    } else {
        event.put(std::move(perThreadStore->getG4Status()));
        event.put(std::move(simsToCheck));
        perThreadStore->putSensitiveDetectorData(simProductGetter);
        perThreadStore->putCutsData(simProductGetter);
    }
    
    
    //perThreadStore->clear();
    
}//end Mu2eG4::produce


// Tell G4 that this run is over.
void Mu2eG4::endRun(art::Run & run, art::ProcessingFrame const& procFrame){

        
        BeamOnEndRun();
        G4cout << "at endRun: numExcludedEvents = " << numExcludedEvents << G4endl;
        
        masterThread->endRun();

}


void Mu2eG4::endJob(art::ProcessingFrame const& procFrame){

    if ( _exportPDTEnd ) exportG4PDT( "End:" );
    physVolHelper_.endRun();
}





/**************************************************************
                    FUNCTION DEFINITIONS
 **************************************************************/


// Do the "end of run" parts of DoEventLoop and BeamOn.
void Mu2eG4::BeamOnEndRun(){

    std::cout << "We are at BeamOnEndRun" <<std::endl;
    //dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->Mu2eG4RunTermination();
  
}//BeamOnEndRun


} // End of namespace mu2e

DEFINE_ART_MODULE(mu2e::Mu2eG4);
