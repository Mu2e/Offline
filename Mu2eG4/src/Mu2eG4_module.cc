// A Producer Module that runs Geant4 and adds its output to the event.
//
// Original author Rob Kutschke
//
// Notes:



// Mu2e includes
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/Mu2eStudyWorld.hh"
#include "Mu2eG4/inc/IMu2eG4Cut.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/exportG4PDT.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eG4/inc/ActionInitialization.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/preG4InitializeTasks.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eG4/inc/generateFieldMap.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"
#include "Mu2eG4/inc/Mu2eG4MultiStageParameters.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/checkConfigRelics.hh"
#include "Mu2eG4/inc/GenEventBroker.hh"
#include "Mu2eG4/inc/EventStash.hh"
#include "Mu2eG4/inc/Mu2eG4MTRunManager.hh"
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined G4VIS_USE_OPENGLQT )
#include "Mu2eG4/inc/Mu2eVisCommands.hh"
#endif

// Data products that will be produced by this module.
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
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"

// Geant4 includes
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined G4VIS_USE_OPENGLQT )
#include "G4VisExecutive.hh"
#endif
#include "G4Run.hh"
#include "G4Timer.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ParticleHPManager.hh"
#include "G4HadronicProcessStore.hh"
#include "G4RunManagerKernel.hh"
#include "G4SDManager.hh"

// C++ includes.
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <utility>

using namespace std;

namespace {
  art::InputTag const invalid_tag{};
}

namespace mu2e {

  class Mu2eG4 : public art::EDProducer {
  public:
    Mu2eG4(fhicl::ParameterSet const& pSet);
    // Accept compiler supplied d'tor

  private:
    void produce(art::Event& e) override;
    void endJob() override;
    void beginRun(art::Run &r) override;
    void endRun(art::Run &) override;
    void beginSubRun(art::SubRun &sr) override;

        fhicl::ParameterSet pset_;

        Mu2eG4ResourceLimits mu2elimits_;
        Mu2eG4TrajectoryControl trajectoryControl_;
        Mu2eG4MultiStageParameters multiStagePars_;

        // The THREE functions that call new G4RunManger functions and break G4's BeamOn() into 3 pieces
        void BeamOnBeginRun( unsigned int runNumber);
        void BeamOnDoOneArtEvent( int eventNumber, G4int, const char* macroFile=0, G4int n_select=-1 );
        void BeamOnEndRun();

        //we need this for MT mode before art-MT is available
        //it fixes the bookkeeping on the art::Ptrs which was messed up by the introduction of the GenParticleStash
        void ReseatPtrsAndMoveDataToArtEvent( art::Event& evt, art::EDProductGetter const* sim_prod_getter, std::unique_ptr<SimParticleCollection> sims_to_store );

        void DoVisualizationFromMacro();

        std::unique_ptr<G4RunManager> _runManager;

        const bool _use_G4MT;
        const G4int _nThreads;

        // Do we issue warnings about multiple runs?
        bool _warnEveryNewRun;

        // Do we want to export the G4 particle data table.
        bool  _exportPDTStart;
        bool  _exportPDTEnd;

        // to be able to make StorePhysicsTable call after the event loop started
        G4VUserPhysicsList* physicsList_;
        std::string storePhysicsTablesDir_;

        ActionInitialization const * _actionInit;

        //these cut objects are used in the master thread to indicate what data product is produced
        //additional thread-local cut objects are owned by ActionInitialization
        std::unique_ptr<IMu2eG4Cut> stackingCuts_;
        std::unique_ptr<IMu2eG4Cut> steppingCuts_;
        std::unique_ptr<IMu2eG4Cut> commonCuts_;

        G4UIsession  *_session;
        G4UImanager  *_UI;
#if     ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined G4VIS_USE_OPENGLQT )
        std::unique_ptr<G4VisManager> _visManager;
#endif

        int _rmvlevel;
        int _tmvlevel;
        int _smvlevel;
        int _checkFieldMap;

        // Names of macro files for visualization.
        string _visMacro;  // init
        string _visGUIMacro; // end of Event GUI

        // Name of a macro file to be used for controling G4 parameters after
        // the initialization phase.
        string _g4Macro;

        art::InputTag _generatorModuleLabel;

        // Helps with indexology related to persisting G4 volume information.
        // string to ptr maps, speed optimization
        // if in MT mode, only allow lookup, don't allow add
        // do a counter that counts how mnay times it was called with an unknown process
        PhysicalVolumeHelper _physVolHelper;

        // handles per-thread objects
        GenEventBroker _genEventBroker;

        // Instance name of the timeVD StepPointMC data product.
        const StepInstanceName _tvdOutputName;
        std::vector<double> timeVDtimes_;

        // Do the G4 initialization that must be done only once per job, not once per run
        void initializeG4( GeometryService& geom, art::Run const& run );

        std::unique_ptr<G4Timer> _timer; // local Mu2e per Geant4 event timer
        // Counters for cumulative time spent processing events by Geant4
        G4double _realElapsed;
        G4double _systemElapsed;
        G4double _userElapsed;

        const bool    _standardMu2eDetector;
        G4ThreeVector _originInWorld;

        std::vector< SensitiveDetectorHelper > SensitiveDetectorHelpers;

        // In sequential mode, there is one SensitiveDetetorHelper object and it is
        // at index 0; in MT mode the SensitiveDetectorHelper object for the master
        // thread is at index _nThreads.
        size_t _masterThreadIndex;

        EventStash _StashForEventData;
        int stashInstanceToStore;

        //used in testing code
        //int event_counter = 0;
        int numExcludedEvents = 0;

  }; // end G4 header

  Mu2eG4::Mu2eG4(fhicl::ParameterSet const& pSet):
    EDProducer{pSet},
    pset_(pSet),
    mu2elimits_(pSet.get<fhicl::ParameterSet>("ResourceLimits")),
    trajectoryControl_(pSet.get<fhicl::ParameterSet>("TrajectoryControl")),
    multiStagePars_(pSet.get<fhicl::ParameterSet>("MultiStageParameters")),

    _use_G4MT(pSet.get<bool>("runinMTMode",false)),
    _nThreads(pSet.get<int>("numberOfThreads",1)),

    _warnEveryNewRun(pSet.get<bool>("debug.warnEveryNewRun",false)),
    _exportPDTStart(pSet.get<bool>("debug.exportPDTStart",false)),
    _exportPDTEnd(pSet.get<bool>("debug.exportPDTEnd",false)),

    storePhysicsTablesDir_(pSet.get<std::string>("debug.storePhysicsTablesDir","")),

    stackingCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4StackingOnlyCut", {}), mu2elimits_)),
    steppingCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4SteppingOnlyCut", {}), mu2elimits_)),
    commonCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4CommonCut", {}), mu2elimits_)),

    _session(nullptr),
    _UI(nullptr),
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined G4VIS_USE_OPENGLQT )
    _visManager(nullptr),
#endif
    // FIXME:  naming of pset parameters
    _rmvlevel(pSet.get<int>("debug.diagLevel",0)),
    _tmvlevel(pSet.get<int>("debug.trackingVerbosityLevel",0)),
    _smvlevel(pSet.get<int>("debug.steppingVerbosityLevel",0)),
    _checkFieldMap(pSet.get<int>("debug.checkFieldMap",0)),
    _visMacro(pSet.get<std::string>("visualization.initMacro")),
    _visGUIMacro(pSet.get<std::string>("visualization.GUIMacro")),
    _g4Macro(pSet.get<std::string>("g4Macro","")),
    _generatorModuleLabel(pSet.get<std::string>("generatorModuleLabel", "")),
    _physVolHelper(),
    //_printPhysicsProcessSummary(pSet.get<bool>("debug.printPhysicsProcessSummary",false)),
    _genEventBroker(_use_G4MT),
    _tvdOutputName(StepInstanceName::timeVD),
    timeVDtimes_(pSet.get<std::vector<double> >("SDConfig.TimeVD.times")),
    _timer(std::make_unique<G4Timer>()),
    _realElapsed(0.),
    _systemElapsed(0.),
    _userElapsed(0.),
    _standardMu2eDetector((art::ServiceHandle<GeometryService>())->isStandardMu2eDetector()),
    _masterThreadIndex(_use_G4MT ? _nThreads : 0),
    _StashForEventData(pSet),
    stashInstanceToStore(-1)
    {
        //get the right type of RunManager depending on whether we are in MT mode or not
        if (_use_G4MT) {
            _runManager.reset(new Mu2eG4MTRunManager());
        }
        else {
            _runManager.reset(new G4RunManager());
        }

        if((_generatorModuleLabel == art::InputTag()) && multiStagePars_.genInputHits().empty()) {
            throw cet::exception("CONFIG")
            << "Error: both generatorModuleLabel and genInputHits are empty - nothing to do!\n";
    }

    //In sequential mode, we need 1 SDHelper and it lives at index 0.
    //
    //In MT, mode we need one SDHelper for each Worker thread, plus one extra for the Master
    //in the ActionInitialization, each worker thread is given one of these SDHs to hold its SD-related data
    //the "0th" worker thread gets the "0th" element of the vector, etc
    //we give the "_nThreads" element to the Master thread through Mu2eG4World to setup the InstanceMap in the ctor of the SDH class
    //we need only one of these SDHs to declare to art the list of products that will be produced
    int nSDHelpersNeeded = ( _use_G4MT ) ? _nThreads+1 : 1;
    SensitiveDetectorHelpers.reserve(nSDHelpersNeeded);

    auto sd_pSet = pSet.get<fhicl::ParameterSet>("SDConfig", fhicl::ParameterSet());;
    for (int i = 0; i != nSDHelpersNeeded; ++i) {
      SensitiveDetectorHelpers.emplace_back(sd_pSet);
    }

    auto& collector = producesCollector();
    SensitiveDetectorHelpers.at(_masterThreadIndex).declareProducts(collector);

    produces<StatusG4>();
    produces<SimParticleCollection>();

    if(!timeVDtimes_.empty()) {
      produces<StepPointMCCollection>(_tvdOutputName.name());
    }

    if(trajectoryControl_.produce()) {
      produces<MCTrajectoryCollection>();
    }

    if(multiStagePars_.multiStage()) {
      produces<SimParticleRemapping>();
    }

    //can we simplify this and directly declare the relevent products
    //rather than contructing these unneccesary object?
    stackingCuts_->declareProducts(collector);
    steppingCuts_->declareProducts(collector);
    commonCuts_->declareProducts(collector);

    // Declare which products this module will read.
    auto const& inputPhysVolTag = multiStagePars_.inputPhysVolumeMultiInfo();
    if (inputPhysVolTag != invalid_tag) {
      consumes<PhysicalVolumeInfoMultiCollection, art::InSubRun>(inputPhysVolTag);
    }
    auto const& inputSimParticlesTag = multiStagePars_.inputSimParticles();
    if (inputSimParticlesTag != invalid_tag) {
      consumes<SimParticleCollection>(inputSimParticlesTag);
    }
    auto const& inputMCTrajectoryTag = multiStagePars_.inputMCTrajectories();
    if (inputMCTrajectoryTag != invalid_tag) {
      consumes<MCTrajectoryCollection>(inputMCTrajectoryTag);
    }
    if (_generatorModuleLabel != invalid_tag) {
      consumes<GenParticleCollection>(_generatorModuleLabel);
    }
    for (auto const& tag : multiStagePars_.genInputHits()) {
      consumes<StepPointMCCollection>(tag);
    }

    produces<PhysicalVolumeInfoMultiCollection,art::InSubRun>();

    // The string "G4Engine" is magic; see the docs for RandomNumberGenerator.
    createEngine( art::ServiceHandle<SeedService>()->getSeed(), "G4Engine");

} // end G4:G4(fhicl::ParameterSet const& pSet);


  // That should really be beginJob().  G4 does not care about run
  // numbers, so we could use a hardcoded 1 for that.  The problem is
  // that Mu2e GeometryService refuses to give information outside of
  // an art run.  That makes sense for alignments and such, but
  // necessitates workarounds for G4 geometry.
void Mu2eG4::beginRun( art::Run &run){

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


    // Tell G4 that we are starting a new run.
    BeamOnBeginRun( run.id().run() );

    // A few more things that only need to be done only once per job,
    // not once per run, but which need to be done after the call to
    // BeamOnReadyToBeginRun.

    //since the cuts used by the individual threads, do we need to do this?
    if ( ncalls == 1 ) {
      stackingCuts_->finishConstruction(_originInWorld);
      steppingCuts_->finishConstruction(_originInWorld);
      commonCuts_->finishConstruction  (_originInWorld);

        //can only be run in single-threaded mode
      if( _checkFieldMap>0 && !(_use_G4MT)) generateFieldMap(_originInWorld,_checkFieldMap);

      if ( _exportPDTStart ) exportG4PDT( "Start:" );//once per job
    }
}


void Mu2eG4::initializeG4( GeometryService& geom, art::Run const& run ){

    //if running in MT mode, set number of threads.
    //need to downcast the ptr to the RunManager, which was defined as a ptr to G4RunManager
    if (_use_G4MT) {
        dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->SetNumberOfThreads(_nThreads);
    }

    if ( _rmvlevel > 0 ) {
      mf::LogInfo logInfo("GEOM");
      logInfo << "Initializing Geant4 for " << run.id()
              << " with verbosity " << _rmvlevel << endl;
      logInfo << " Configured simParticleNumberOffset = "<< multiStagePars_.simParticleNumberOffset() << endl;
    }


    // Create user actions and register them with G4.
    G4VUserDetectorConstruction* allMu2e;

    // as mentioned above, we give the last element to the Master thread to setup the InstanceMap in the ctor of the SDH class
    
    if (_standardMu2eDetector) {
      geom.addWorldG4(*GeomHandle<Mu2eHall>());

      allMu2e = 
     	(new WorldMaker<Mu2eWorld>(std::make_unique<Mu2eWorld>(pset_, &(SensitiveDetectorHelpers.at(_masterThreadIndex))  ),
     				   std::make_unique<ConstructMaterials>(pset_)) );
    
      _originInWorld = (GeomHandle<WorldG4>())->mu2eOriginInWorld();
    }
    else {
      allMu2e =
	(new WorldMaker<Mu2eStudyWorld>(std::make_unique<Mu2eStudyWorld>(pset_, &(SensitiveDetectorHelpers.at(_masterThreadIndex)) ),
					std::make_unique<ConstructMaterials>(pset_)) );

      // non-Mu2e detector: the system origin os set to (0.,0.,0.); do not use geometry service for that
      _originInWorld = G4ThreeVector(0.0,0.0,0.0);
    }


    // in the non Mu2e detector we are working in the system with the
    // origin set to 0.,0.,0. and do not use geometry service for that
    //    originInWorld = (!_standardMu2eDetector) ? G4ThreeVector(0.0,0.0,0.0) : (GeomHandle<WorldG4>())->mu2eOriginInWorld();


    preG4InitializeTasks(pset_);

    _runManager->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(allMu2e);

    physicsList_ = physicsListDecider(pset_);
    physicsList_->SetVerboseLevel(_rmvlevel);

    G4ParticleHPManager::GetInstance()->SetVerboseLevel(_rmvlevel);

    G4HadronicProcessStore::Instance()->SetVerbose(_rmvlevel);

    _runManager->SetUserInitialization(physicsList_);


    //this is where the UserActions are instantiated
    ActionInitialization* actioninit = new ActionInitialization(pset_,
                                                                SensitiveDetectorHelpers,
                                                                &_genEventBroker, &_physVolHelper,
                                                                _use_G4MT, _nThreads, _originInWorld,
                                                                mu2elimits_,
                                                                multiStagePars_.simParticleNumberOffset()
                                                                );

    //in MT mode, this is where BuildForMaster is called for master thread
    // in sequential mode, this is where Build() is called for main thread
    _runManager->SetUserInitialization(actioninit);

    _UI = G4UImanager::GetUIpointer();

    // Any final G4 interactive commands ...
    if ( !_g4Macro.empty() ) {
        G4String command("/control/execute ");
        ConfigFileLookupPolicy path;
        command += path(_g4Macro);
        _UI->ApplyCommand(command);

    }

    if ( _rmvlevel > 0 ) {

      //  GetRunManagerType() returns an enum named RMType to indicate
      //  what kind of RunManager it is. RMType is defined as {
      //  sequentialRM, masterRM, workerRM }

       G4cout << __func__
              << " Before Initialize(), G4RunManagerType: "
              << _runManager->GetRunManagerType() << G4endl;
    }
    // Initialize G4 for this run.
    _runManager->Initialize();

    if ( _rmvlevel > 0 ) {
       G4cout << __func__
              << " After Initialize(), G4RunManagerType:  "
              << _runManager->GetRunManagerType() << G4endl;
    }

#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined G4VIS_USE_OPENGLQT )
    // Setup the graphics if requested.
    if ( !_visMacro.empty() ) {

        _visManager = std::unique_ptr<G4VisManager>(new G4VisExecutive());
        _visManager->Initialize();

        ConfigFileLookupPolicy visPath;

        G4String command("/control/execute ");
        command += visPath(_visMacro);

        _UI->ApplyCommand( command );

    }
#endif

} // end G4::initializeG4


void Mu2eG4::beginSubRun(art::SubRun& sr)
{
    using Collection_t = PhysicalVolumeInfoMultiCollection;
    auto mvi = std::make_unique<Collection_t>();

    if(multiStagePars_.inputPhysVolumeMultiInfo() != invalid_tag) {
      // Copy over data from the previous simulation stages
      auto const& ih = sr.getValidHandle<Collection_t>(multiStagePars_.inputPhysVolumeMultiInfo());
      mvi->reserve(1 + ih->size());
      mvi->insert(mvi->begin(), ih->cbegin(), ih->cend());

    }

    // Append info for the current stage
    mvi->emplace_back(multiStagePars_.simParticleNumberOffset(), _physVolHelper.persistentSingleStageInfo());

    sr.put(std::move(mvi));
}


// Create one G4 event and copy its output to the art::event.
void Mu2eG4::produce(art::Event& event) {

    if (_use_G4MT && event.id().event() == 1) std::cout << "\n*-*-*-*-*-*- You are running "
       << dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->GetNumberOfThreads() << " threads. -*-*-*-*-*-*\n" << std::endl;

    //confirm that IF we are running in MT mode we do not have inputs from previous simulation stages
    //otherwsie, throw an exception
    if (_use_G4MT) {

        if (   multiStagePars_.inputSimParticles() != art::InputTag()
            || multiStagePars_.inputMCTrajectories() != art::InputTag()
            || !(multiStagePars_.genInputHits().empty()) ) {
            throw cet::exception("CONFIG")
            << "Error: You are trying to run in MT mode with input from previous stages.  This is an invalid configuration!\n";
        }
    }

    // ProductID and ProductGetter for the SimParticleCollection.
    art::ProductID simPartId(event.getProductID<SimParticleCollection>());
    art::EDProductGetter const* simProductGetter = event.productGetter(simPartId);

    //stash is empty, we need to simulate events
    if (_StashForEventData.getStashSize() == 0)
    {

        stashInstanceToStore = 0;

        //******** these are per art::event quantities ********
        // StepPointMCCollection of input hits from the previous simulation stage
        HitHandles genInputHits;
        for(const auto& i : multiStagePars_.genInputHits()) {
            genInputHits.emplace_back(event.getValidHandle<StepPointMCCollection>(i));
        }

        _genEventBroker.loadEvent(genInputHits, simPartId, &event, _generatorModuleLabel, &_StashForEventData, simProductGetter);
        //getStashSize() can only be called after loadEvent is called

        if (_use_G4MT)//in MT mode, stash size is given by the size of input GenParticleCollection
        {
            if (_genEventBroker.getStashSize() == 0) {
                throw cet::exception("CONFIG")
                << "Error: You are trying to run in MT mode with a stash size of '0'.  This is an invalid configuration!\n";
            }

            _StashForEventData.initializeStash(_genEventBroker.getStashSize());
        }
        else//in sequential mode, the stash size is 1
        {
            _StashForEventData.initializeStash(1);
        }

        // Run G4 for this event and access the completed event.
        BeamOnDoOneArtEvent( event.id().event(), _genEventBroker.getStashSize() );

        _genEventBroker.setEventPtrToZero();


    }//end if stash is empty, simulate events


    //testing stuff ********************************
    //    std::cout << "in produce, printing the Stash Sim Particle info " << std::endl;
    //    _StashForEventData.printInfo(stashInstanceToStore);


    std::unique_ptr<SimParticleCollection> simsToCheck = std::move(_StashForEventData.getSimPartCollection(stashInstanceToStore));

    if (simsToCheck == nullptr) {
        numExcludedEvents++;
    } else {

        if (_use_G4MT)
        {
            //now move the data into the art::Event IFF the event has passed the StepPointMomentumFilter inside Mu2eG4EventAction
            //if event has not passed, the ptr to the SimParticleCollection will be null
            event.put(std::move(_StashForEventData.getG4Status(stashInstanceToStore)));
            ReseatPtrsAndMoveDataToArtEvent(event, simProductGetter, std::move(simsToCheck));
            _StashForEventData.putSensitiveDetectorData(stashInstanceToStore, event, simProductGetter);
            _StashForEventData.putCutsData(stashInstanceToStore, event, simProductGetter);

        }
        else//in sequential mode, there isn't any need to reseat the ptrs, since there is no GenParticleCollections object
        {

            event.put(std::move(_StashForEventData.getG4Status(stashInstanceToStore)));
            event.put(std::move(std::move(simsToCheck)));

            if(!timeVDtimes_.empty()) {
                event.put(std::move(_StashForEventData.getTVDHits(stashInstanceToStore)),_StashForEventData.getTVDName(stashInstanceToStore));
            }

            if(trajectoryControl_.produce()) {
                event.put(std::move(_StashForEventData.getMCTrajCollection(stashInstanceToStore)));
            }

            if(multiStagePars_.multiStage()) {
                event.put(std::move(_StashForEventData.getSimParticleRemap(stashInstanceToStore)));
            }

            if(SensitiveDetectorHelpers[_masterThreadIndex].extMonPixelsEnabled()) {
                event.put(std::move(_StashForEventData.getExtMonFNALSimHitCollection(stashInstanceToStore)));
            }

            _StashForEventData.putSensitiveDetectorData(stashInstanceToStore, event, simProductGetter);
            _StashForEventData.putCutsData(stashInstanceToStore, event, simProductGetter);
        }//sequential
    }//simsToCheck !=nullptr

    //increment the instance of the EventStash to store
    stashInstanceToStore++;

    if (stashInstanceToStore == _StashForEventData.getStashSize()) {
        _StashForEventData.clearStash();
    }


}//end Mu2eG4::produce


// Tell G4 that this run is over.
void Mu2eG4::endRun(art::Run & run){

        BeamOnEndRun();

        G4cout << "at endRun: numExcludedEvents = " << numExcludedEvents << G4endl;

}


void Mu2eG4::endJob(){

    if ( _exportPDTEnd ) exportG4PDT( "End:" );

    // Yes, these are named endRun, but they are really endJob actions.
    _physVolHelper.endRun();


    //this PhysicsProcessInfo has been made a per-thread item
   // G4AutoLock PIlock(&processInfoMutex);
    //if ( _printPhysicsProcessSummary ){ _processInfo.endRun(); }

    //_processInfo.endRun();

    //PIlock.unlock();
}





/**************************************************************
                    FUNCTION DEFINITIONS
 **************************************************************/

// Do the "begin run" parts of BeamOn.
void Mu2eG4::BeamOnBeginRun( unsigned int runNumber){

  if ( _rmvlevel > 0 ) {
    G4cout << __func__ << ": runNumber: " << runNumber << G4endl;
  }


        _runManager->SetRunIDCounter(runNumber);

        bool cond = _runManager->ConfirmBeamOnCondition();
        if(!cond){
            // throw here
            return;
        }

        _realElapsed   = 0.;
        _systemElapsed = 0.;
        _userElapsed   = 0.;

        //NOTE: the data member G4RunManager::numberOfEventToBeProcessed would have been set from within the G4 call to BeamOn()
        //BEFORE the call the ConstructScoringWorlds() with this call:
        //_runManager->SetNumberOfEventsToBeProcessed(numberOfEventsToBeProcessed);
        //Since we have deconstructed the BeamOn() function, we must set this data member ourselves.
        //In the pre-MT version of this code, this call was made here, immediately before the call _runManager->ConstructScoringWorlds();
        //It has been moved, because at 'BeginRun' time, the value of this variable cannot be known without using a FHiCL parameter.
        //It was decided not to define it though a FHiCL parameter, but rather within the code.
        //
        //In SEQUENTIAL mode, this variable doesn't seem to serve any purpose other than to be used in the call
        //currentRun->SetNumberOfEventToBeProcessed(numberOfEventToBeProcessed); in G4RunManager::RunInitialization().
        //Furthermore, at least up through G4 v10.4.02 the data member G4Run::numberOfEventToBeProcessed is only informative.
        //It is not needed anywhere in the G4 code.  So, it is OK to not set this variable at this point.
        //
        //For MT mode, G4MTRunManager::numberOfEventToBeProcessed is used to distribute events to the threads, among other things,
        //and must have the value of EITHER the size of the GenParticle Stash when using single-threaded art,
        //OR the total number of events to be simulated when using multithreaded art.
        //It must be set before the call to G4MTRunManager::InitializeEventLoop.
        //We can get the size of the GenParticle Stash from WITHIN the code with the call to _genEventBroker.getStashSize(),
        //but this cannot be called until after _genEventBroker.loadEvent() is called.
        //So, we get this variable from the GenParticle Stash after the call to _genEventBroker.loadEvent().  We initialize the
        //EventStash to the correct size using it, and we pass it into G4 in the call BeamOnDoOneArtEvent().

        _runManager->ConstructScoringWorlds();
        _runManager->RunInitialization();

}//BeamOnBeginRun


// Do the "per event" part of DoEventLoop.
void Mu2eG4::BeamOnDoOneArtEvent( int eventNumber, G4int num_events, const char* macroFile, G4int n_select){

  if ( _rmvlevel > 1 ) {
    G4cout << __func__ << ": eventNumber: " << eventNumber << ", num_events: " << num_events << G4endl;
  }



        if (_use_G4MT)//MT mode
        {
            //this is where the events are actually processed
            //num_events is # of G4 events processed per art event

            //NOTE: G4MTRunManager::WaitForEndEventLoopWorkers() is a protected function in G4 code, so we MUST have our own MTRunManager
            //in order to access this function

            dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->SetNumberOfEventsToBeProcessed(num_events);
            dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->InitializeEventLoop(num_events,macroFile,n_select);
            dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->Mu2eG4WaitForEndEventLoopWorkers();
            //dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->Mu2eG4TerminateWorkers();//DID NOT WORK BEFORE

        }
        else//sequential mode
        {
            _runManager->InitializeEventLoop(num_events,macroFile,n_select);

            _timer->Start();
            _runManager->ProcessOneEvent(eventNumber);
            _timer->Stop();

            // Accumulate time spent in G4 for all events in this run.
            _realElapsed   += _timer->GetRealElapsed();
            _systemElapsed += _timer->GetSystemElapsed();
            _userElapsed   += _timer->GetUserElapsed();


            // Pause to see graphics.
            if ( !_visMacro.empty() ){
                DoVisualizationFromMacro();
            }
            _runManager->TerminateOneEvent();
        }//end if

}//BeamOnDoOneArtEvent


// Do the "end of run" parts of DoEventLoop and BeamOn.
void Mu2eG4::BeamOnEndRun(){

    if (storePhysicsTablesDir_!="") {
      if ( _rmvlevel > 0 ) {
        G4cout << __func__ << " Will write out physics tables to "
               << storePhysicsTablesDir_
               << G4endl;
      }
      physicsList_->StorePhysicsTable(storePhysicsTablesDir_);
    }

    if (_use_G4MT)//MT mode
    {
        dynamic_cast<Mu2eG4MTRunManager*>(_runManager.get())->Mu2eG4RunTermination();
    }
    else//sequential mode
    {
        _runManager->TerminateEventLoop();
        _runManager->RunTermination();

        G4cout << "  Event processing inside ProcessOneEvent time summary" << G4endl;
        G4cout << "  User="  << _userElapsed
        << "s Real="  << _realElapsed
        << "s Sys="   << _systemElapsed
        << "s" << G4endl;
    }
}//BeamOnEndRun

void Mu2eG4::ReseatPtrsAndMoveDataToArtEvent(art::Event& evt, art::EDProductGetter const* sim_prod_getter, std::unique_ptr<SimParticleCollection> sims_to_store){

    art::Handle<GenParticleCollection> gensHandle;
    if(!(_generatorModuleLabel == art::InputTag())) {
        evt.getByLabel(_generatorModuleLabel, gensHandle);
    }

    //put the SimParticleCollection into the event
    //std::unique_ptr<SimParticleCollection> tempSims = std::move(_StashForEventData.getSimPartCollection(stashInstanceToStore));

    for ( SimParticleCollection::iterator i=sims_to_store->begin(); i!=sims_to_store->end(); ++i )
    {
        SimParticle& sim = i->second;

        if ( _use_G4MT && sim.isPrimary() && gensHandle.isValid() ){
            art::Ptr<GenParticle> reseat(gensHandle, sim.genParticle().key());
            sim.genParticle() = reseat;
        }

        sim.parent() = art::Ptr<SimParticle>(sim.parent().id(),
                                             sim.parent().key(),
                                             sim_prod_getter );

        std::vector<art::Ptr<SimParticle> > const& daughters = sim.daughters();

        if ( !daughters.empty() ) {
            std::vector<art::Ptr<SimParticle> > newDaughters;
            newDaughters.reserve(daughters.size());

            for ( size_t i=0; i != daughters.size(); ++i){
                art::Ptr<SimParticle> const& dau = art::Ptr<SimParticle>(daughters[i].id(), daughters[i].key(),
                                                                         sim_prod_getter );
                newDaughters.push_back( dau );
            }
            sim.setDaughterPtrs( newDaughters );
        }

    }//for (SimParticleCollection::iterator...

    evt.put(std::move(sims_to_store));


    if(!timeVDtimes_.empty()) {
        std::unique_ptr<StepPointMCCollection> tempTVD = std::move(_StashForEventData.getTVDHits(stashInstanceToStore));

        for ( StepPointMCCollection::iterator i=tempTVD->begin(); i!=tempTVD->end(); ++i ){
            StepPointMC& step = *i;

            if ( step.simParticle().isNonnull() ){
                step.simParticle() = art::Ptr<SimParticle>(step.simParticle().id(),
                                                           step.simParticle().key(),
                                                           sim_prod_getter );
            }
        }
        evt.put(std::move(tempTVD),_StashForEventData.getTVDName(stashInstanceToStore));
    }// if !timeVDtimes_.empty()


    if(trajectoryControl_.produce()) {
        //get the MCTrajCollection from the Stash and create a new one to put stuff into
        std::unique_ptr<MCTrajectoryCollection> tempTrajs = std::move(_StashForEventData.getMCTrajCollection(stashInstanceToStore));
        std::unique_ptr<MCTrajectoryCollection> outTrajectory(new MCTrajectoryCollection());

        for ( MCTrajectoryCollection::iterator i=tempTrajs->begin(); i!=tempTrajs->end(); ++i ){
            art::Ptr<SimParticle> newParticle(i->second.sim().id(), i->second.sim().key(), sim_prod_getter );
            (*outTrajectory)[newParticle] = i->second;
            (*outTrajectory)[newParticle].sim() = newParticle;
        }
        evt.put(std::move(outTrajectory));
    }// if trajectoryControl


    //I am including this here for symmetry.  These are not produced in MT mode.
    //DO I NEED TO RESEAT THE SimParticles in the Remap?  No, becasue these are only produced in sequential mode
    //where the stash operates on a one-event-in, one-event-out basis
    if(multiStagePars_.multiStage()) {
        evt.put(std::move(_StashForEventData.getSimParticleRemap(stashInstanceToStore)));
    }

    // Fixme: does this work in MT mode?  If not, does it need to?
    if(SensitiveDetectorHelpers.at(_masterThreadIndex).extMonPixelsEnabled()) {
        std::unique_ptr<ExtMonFNALSimHitCollection> tempExtMonHits = std::move(_StashForEventData.getExtMonFNALSimHitCollection(stashInstanceToStore));

        for ( ExtMonFNALSimHitCollection::iterator i=tempExtMonHits->begin(); i!=tempExtMonHits->end(); ++i ){
            ExtMonFNALSimHit& hit = *i;

            if ( hit.simParticle().isNonnull() ){
                hit.simParticle() = art::Ptr<SimParticle>(hit.simParticle().id(), hit.simParticle().key(), sim_prod_getter );
            }
        }
        evt.put(std::move(tempExtMonHits));
    }//if extMonPixelsEnabled

}//ReseatPtrsAndMoveDataToArtEvent

void Mu2eG4::DoVisualizationFromMacro(){

    // Prompt to continue and wait for reply.
    cout << "Enter a character to go to the next event" << endl;
    cout << "q quits, s enters G4 interactive session, g enters a GUI session (if available)"
    << endl;
    cout << "Once in G4 interactive session, to quit it type \"exit\" or use File menu"
    << endl;

    string userinput;
    cin >> userinput;
    G4cout << userinput << G4endl;

    // Check if user is requesting an early termination of the event loop.
    if ( !userinput.empty() ){
        // Check only the first character; >> skips whitespace by default
        char c = tolower( userinput[0] );
        if ( c == 'q' ){
            throw cet::exception("CONTROL")
            << "Early end of event loop requested inside G4, \n";
        } else if ( c == 's' || c == 'g' || c == 'v' ){
            // v is for backward compatibility
            G4int argc=1;
            // Cast away const-ness; required by the G4 interface ...
            char* dummy = (char *)"dummy";
            char** argv = &dummy;
            G4UIExecutive* UIE = ( c == 's' || c == 'v' ) ?
            new G4UIExecutive(argc, argv,"tcsh") :
            new G4UIExecutive(argc, argv);

#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined G4VIS_USE_OPENGLQT )

            if (UIE->IsGUI()) {

                // we add a command here and initialize it
                // (/vis/sceneHandler has to exist prior to this)
                auto* drEv = new Mu2eVisCommandSceneHandlerDrawEvent();
                _visManager->RegisterMessenger(drEv); // assumes ownership;
                // drEv->SetVisManager(_visManager.get());
                // vis manager pointer is static member of the drEv base
                // class so the above is not needed

                if ( !_visGUIMacro.empty() ){
                    G4String command("/control/execute ");
                    ConfigFileLookupPolicy visPath;
                    command += visPath(_visGUIMacro);
                    _UI->ApplyCommand( command );

                    cout << "In GUI interactive session use the \"Start Here\" menu "
                    << "followed by the Viewer commands or redisplaying event"
                    << endl;

                } else {
                    cout << __func__ << " WARNING: visGUIMacro empty, may need to be defined in fcl" << endl;
                }

            } // end UIE->IsGUI()
#endif
            UIE->SessionStart();
            delete UIE;

            //If current scene is scene-0 and if scene-handler-0 has viewer-0 we
            //will select it if not current to deal with a case which may occur
            //e.g. in a simultaneous use of OGL & Qt

            // basically _UI->ApplyCommand("/vis/viewer/select viewer-0"); // to have tracks drawn

#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined  G4VIS_USE_OPENGLQT )
            G4String viewerToLookFor("viewer-0");
            G4VViewer* pViewer = _visManager->GetViewer(viewerToLookFor);
            if (pViewer) {
                if (pViewer != _visManager->GetCurrentViewer()) {
                    _visManager->SetCurrentViewer(pViewer);
                }
            }
            // G4VGraphicsSystem* gsys = _visManager->GetCurrentGraphicsSystem();
            // if (gsys) {
            //   cout << __func__ << " current GraphicsSystem Name " << gsys->GetName() <<  endl;
            // }
#endif
        } // end c == 'q'

    } // end !userinput.empty()

}//DoVisualizationFromMacro


} // End of namespace mu2e

DEFINE_ART_MODULE(mu2e::Mu2eG4);
