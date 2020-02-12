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
#include "Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Mu2eG4/inc/Mu2eG4Config.hh"
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined G4VIS_USE_OPENGLQT )
#include "Mu2eG4/inc/Mu2eVisCommands.hh"
#endif

// Data products that will be produced by this module.
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"

// From art and its tool chain.
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
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
#include "G4RunManager.hh"
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

    using Parameters = art::EDProducer::Table<Mu2eG4Config::Top>;
    explicit Mu2eG4(const Parameters& pars);

  private:
    void produce(art::Event& e) override;
    void endJob() override;
    void beginRun(art::Run &r) override;
    void endRun(art::Run &) override;
    void beginSubRun(art::SubRun &sr) override;

    Mu2eG4Config::Top conf_;

    Mu2eG4ResourceLimits mu2elimits_;
    Mu2eG4TrajectoryControl trajectoryControl_;
    Mu2eG4MultiStageParameters multiStagePars_;

    // The THREE functions that call new G4RunManger functions and break G4's BeamOn() into 3 pieces
    void BeamOnBeginRun( unsigned int runNumber, const char* macroFile=0, G4int n_select=-1 );
    void BeamOnDoOneArtEvent( int eventNumber );
    void BeamOnEndRun();

    void DoVisualizationFromMacro();

    std::unique_ptr<G4RunManager> _runManager;

    // Do we issue warnings about multiple runs?
    bool _warnEveryNewRun;

    // Do we want to export the G4 particle data table.
    bool  _exportPDTStart;
    bool  _exportPDTEnd;

    // to be able to make StorePhysicsTable call after the event loop started
    G4VUserPhysicsList* physicsList_;
    std::string storePhysicsTablesDir_;

    //these cut objects are used to indicate what data product is produced
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
    // do a counter that counts how mnay times it was called with an unknown process
    PhysicalVolumeHelper _physVolHelper;

    // Instance name of the timeVD StepPointMC data product.
    const StepInstanceName _tvdOutputName;
    bool timeVD_enabled_;

    // Do the G4 initialization that must be done only once per job, not once per run
    void initializeG4( GeometryService& geom, art::Run const& run );

    std::unique_ptr<G4Timer> _timer; // local Mu2e per Geant4 event timer
    // Counters for cumulative time spent processing events by Geant4
    G4double _realElapsed;
    G4double _systemElapsed;
    G4double _userElapsed;

    const bool    _standardMu2eDetector;
    G4ThreeVector _originInWorld;

    SensitiveDetectorHelper _sensitiveDetectorHelper;

    Mu2eG4PerThreadStorage perThreadStore;

    int numExcludedEvents = 0;

  }; // end G4 header

  Mu2eG4::Mu2eG4(const Parameters& pars):
    EDProducer{pars},
    conf_(pars()),
    mu2elimits_(pars().ResourceLimits()),
    trajectoryControl_(pars().TrajectoryControl()),
    multiStagePars_(pars()),
    _runManager(std::make_unique<G4RunManager>()),
    _warnEveryNewRun(pars().debug().warnEveryNewRun()),
    _exportPDTStart(pars().debug().exportPDTStart()),
    _exportPDTEnd(pars().debug().exportPDTEnd()),

    storePhysicsTablesDir_(pars().debug().storePhysicsTablesDir()),

    stackingCuts_(createMu2eG4Cuts(pars().Mu2eG4StackingOnlyCut.get<fhicl::ParameterSet>(), mu2elimits_)),
    steppingCuts_(createMu2eG4Cuts(pars().Mu2eG4SteppingOnlyCut.get<fhicl::ParameterSet>(), mu2elimits_)),
    commonCuts_(createMu2eG4Cuts(pars().Mu2eG4CommonCut.get<fhicl::ParameterSet>(), mu2elimits_)),

    _session(nullptr),
    _UI(nullptr),
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined G4VIS_USE_OPENGLQT )
    _visManager(nullptr),
#endif
    _rmvlevel(pars().debug().diagLevel()),
    _checkFieldMap(pars().debug().checkFieldMap()),
    _visMacro(pars().visualization().initMacro()),
    _visGUIMacro(pars().visualization().GUIMacro()),
    _g4Macro(pars().g4Macro()),
    _generatorModuleLabel(pars().generatorModuleLabel()),
    _physVolHelper(),
    _tvdOutputName(StepInstanceName::timeVD),
    timeVD_enabled_(pars().SDConfig().TimeVD().enabled()),
    _timer(std::make_unique<G4Timer>()),
    _realElapsed(0.),
    _systemElapsed(0.),
    _userElapsed(0.),
    _standardMu2eDetector((art::ServiceHandle<GeometryService>())->isStandardMu2eDetector()),
    _sensitiveDetectorHelper(pars().SDConfig()),
    perThreadStore()
    {

      if((_generatorModuleLabel == art::InputTag()) && multiStagePars_.genInputHits().empty()) {
        throw cet::exception("CONFIG")
          << "Error: both generatorModuleLabel and genInputHits are empty - nothing to do!\n";
      }

      auto& collector = producesCollector();
      _sensitiveDetectorHelper.declareProducts(collector);

      produces<StatusG4>();
      produces<SimParticleCollection>();

      if(timeVD_enabled_) {
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

    } // end Mu2eG4 constructor


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

    if ( ncalls == 1 ) {
      if( _checkFieldMap>0 ) generateFieldMap(_originInWorld,_checkFieldMap);
      if ( _exportPDTStart ) exportG4PDT( "Start:" );//once per job
    }
  }


  void Mu2eG4::initializeG4( GeometryService& geom, art::Run const& run ){

    if ( _rmvlevel > 0 ) {
      mf::LogInfo logInfo("GEOM");
      logInfo << "Initializing Geant4 for " << run.id()
              << " with verbosity " << _rmvlevel << endl;
      logInfo << " Configured simParticleNumberOffset = "<< multiStagePars_.simParticleNumberOffset() << endl;
    }


    // Create user actions and register them with G4.
    G4VUserDetectorConstruction* allMu2e = nullptr;

    // as mentioned above, we give the last element to the Master thread to setup the InstanceMap in the ctor of the SDH class

    if (_standardMu2eDetector) {
      geom.addWorldG4(*GeomHandle<Mu2eHall>());

      allMu2e =
        (new WorldMaker<Mu2eWorld>(std::make_unique<Mu2eWorld>(conf_, &(_sensitiveDetectorHelper)  ),
                                   std::make_unique<ConstructMaterials>(conf_.debug())) );

      _originInWorld = (GeomHandle<WorldG4>())->mu2eOriginInWorld();
    }
    else {
      allMu2e =
        (new WorldMaker<Mu2eStudyWorld>(std::make_unique<Mu2eStudyWorld>(conf_, &(_sensitiveDetectorHelper) ),
                                        std::make_unique<ConstructMaterials>(conf_.debug())) );

      // non-Mu2e detector: the system origin os set to (0.,0.,0.); do not use geometry service for that
      _originInWorld = G4ThreeVector(0.0,0.0,0.0);
    }

    preG4InitializeTasks(conf_.physics(), conf_.debug());

    _runManager->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(allMu2e);

    physicsList_ = physicsListDecider(conf_.physics(), conf_.debug());
    physicsList_->SetVerboseLevel(_rmvlevel);

    G4ParticleHPManager::GetInstance()->SetVerboseLevel(_rmvlevel);

    G4HadronicProcessStore::Instance()->SetVerbose(_rmvlevel);

    _runManager->SetUserInitialization(physicsList_);


    //this is where the UserActions are instantiated
    ActionInitialization* actioninit = new ActionInitialization(conf_,
                                                                &_sensitiveDetectorHelper,
                                                                &perThreadStore,
                                                                &_physVolHelper,
                                                                _originInWorld,
                                                                multiStagePars_.simParticleNumberOffset()
                                                                );

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

    art::Handle<GenParticleCollection> gensHandle;
    if(!(_generatorModuleLabel == art::InputTag())) {
      event.getByLabel(_generatorModuleLabel, gensHandle);
    }

    // StepPointMCCollection of input hits from the previous simulation stage
    HitHandles genInputHits;
    for(const auto& i : multiStagePars_.genInputHits()) {
      genInputHits.emplace_back(event.getValidHandle<StepPointMCCollection>(i));
    }

    // ProductID and ProductGetter for the SimParticleCollection.
    art::ProductID simPartId(event.getProductID<SimParticleCollection>());
    art::EDProductGetter const* simProductGetter = event.productGetter(simPartId);

    SimParticleHelper spHelper(multiStagePars_.simParticleNumberOffset(), simPartId, &event, simProductGetter);
    SimParticlePrimaryHelper parentHelper(&event, simPartId, gensHandle, simProductGetter);

    perThreadStore.initializeEventInfo(&event, &spHelper, &parentHelper, &genInputHits, _generatorModuleLabel);

    // Run G4 for this event and access the completed event.
    BeamOnDoOneArtEvent( event.id().event() );


    /////////////////////////////////////////////////////////////////////////////////////
    std::unique_ptr<SimParticleCollection> simsToCheck = perThreadStore.getSimPartCollection();

    if (simsToCheck == nullptr) {
      numExcludedEvents++;
    }
    else {

      event.put(std::move(perThreadStore.getG4Status()));
      event.put(std::move(simsToCheck));
      perThreadStore.putSensitiveDetectorData(simProductGetter);
      perThreadStore.putCutsData(simProductGetter);

      if(timeVD_enabled_) {
        event.put(std::move(perThreadStore.getTVDHits()),perThreadStore.getTVDName());
      }

      if(trajectoryControl_.produce()) {
        event.put(std::move(perThreadStore.getMCTrajCollection()));
      }

      if(multiStagePars_.multiStage()) {
        event.put(std::move(perThreadStore.getSimParticleRemap()));
      }

      if(_sensitiveDetectorHelper.extMonPixelsEnabled()) {
        event.put(std::move(perThreadStore.getExtMonFNALSimHitCollection()));
      }

    }//simsToCheck !=nullptr

    _runManager->TerminateOneEvent();
    perThreadStore.clearData();
  }//end Mu2eG4::produce


  // Tell G4 that this run is over.
  void Mu2eG4::endRun(art::Run & run){

    BeamOnEndRun();
    G4cout << "at endRun: numExcludedEvents = " << numExcludedEvents << G4endl;
  }


  void Mu2eG4::endJob(){

    if ( _exportPDTEnd ) exportG4PDT( "End:" );
    _physVolHelper.endRun();
  }





  /**************************************************************
                    FUNCTION DEFINITIONS
  **************************************************************/

  // Do the "begin run" parts of BeamOn.
  void Mu2eG4::BeamOnBeginRun( unsigned int runNumber, const char* macroFile, G4int n_select){

    if ( _rmvlevel > 0 ) {
      G4cout << __func__ << ": runNumber: " << runNumber << G4endl;
    }

    _runManager->SetRunIDCounter(runNumber);

    bool cond = _runManager->ConfirmBeamOnCondition();
    if(!cond){
      throw cet::exception("G4INIT")
        << "ConfirmBeamOn failed\n";
    }

    _realElapsed   = 0.;
    _systemElapsed = 0.;
    _userElapsed   = 0.;


    G4int numberOfEventsToBeProcessed = std::numeric_limits<int>::max(); // largest int for now

    _runManager->SetNumberOfEventsToBeProcessed(numberOfEventsToBeProcessed);// this would have been set by BeamOn
    _runManager->ConstructScoringWorlds();
    _runManager->RunInitialization();

    _runManager->InitializeEventLoop(numberOfEventsToBeProcessed,macroFile,n_select);


  }//BeamOnBeginRun


  // Do the "per event" part of DoEventLoop.
  void Mu2eG4::BeamOnDoOneArtEvent( int eventNumber ){

    if ( _rmvlevel > 1 ) {
      G4cout << __func__ << ": eventNumber: " << eventNumber << G4endl;
    }

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

    _runManager->TerminateEventLoop();
    _runManager->RunTermination();

    G4cout << "  Event processing inside ProcessOneEvent time summary" << G4endl;
    G4cout << "  User="  << _userElapsed
           << "s Real="  << _realElapsed
           << "s Sys="   << _systemElapsed
           << "s" << G4endl;

  }//BeamOnEndRun


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
