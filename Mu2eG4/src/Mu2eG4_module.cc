// A Producer Module that runs Geant4 and adds its output to the event.
//
// Original author Rob Kutschke
//
// Notes:



// Mu2e includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Offline/Mu2eG4/inc/WorldMaker.hh"
#include "Offline/Mu2eG4/inc/Mu2eWorld.hh"
#include "Offline/Mu2eG4/inc/Mu2eStudyWorld.hh"
#include "Offline/Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Offline/Mu2eG4/inc/exportG4PDT.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ActionInitialization.hh"
#include "Offline/Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Offline/Mu2eG4/inc/physicsListDecider.hh"
#include "Offline/Mu2eG4/inc/preG4InitializeTasks.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4SensitiveDetector.hh"
#include "Offline/Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/Mu2eG4/inc/generateFieldMap.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ResourceLimits.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4TrajectoryControl.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4Inputs.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/checkConfigRelics.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4IOConfigHelper.hh"
#include "Offline/Mu2eG4/inc/validGeometryOrThrow.hh"
#include "Offline/Mu2eG4/inc/writePhysicalVolumes.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4ScoringManager.hh"
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined G4VIS_USE_OPENGLQT )
#include "Offline/Mu2eG4/inc/Mu2eG4VisCommands.hh"
#endif

// Data products that will be produced by this module.
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Offline/MCDataProducts/inc/StatusG4.hh"
#include "Offline/MCDataProducts/inc/StepInstanceName.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/SimParticleRemapping.hh"

// From art and its tool chain.
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"

// Geant4 includes
#include "Geant4/G4UIExecutive.hh"
#include "Geant4/G4UImanager.hh"
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined G4VIS_USE_OPENGLQT )
#include "Geant4/G4VisExecutive.hh"
#endif
#include "Geant4/G4Run.hh"
#include "Geant4/G4Timer.hh"
#include "Geant4/G4VUserPhysicsList.hh"
#if G4VERSION>4106
#include "Geant4/G4HadronicParameters.hh"
#else
#include "Geant4/G4ParticleHPManager.hh"
#include "Geant4/G4HadronicProcessStore.hh"
#endif
#include "Geant4/G4RunManagerKernel.hh"
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4ScoringManager.hh"

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
    void endSubRun(art::SubRun &sr) override;

    Mu2eG4Config::Top conf_;

    Mu2eG4ResourceLimits mu2elimits_;
    Mu2eG4TrajectoryControl trajectoryControl_;
    Mu2eG4Inputs multiStagePars_;

    unsigned simStage_;

    // The THREE functions that call new G4RunManger functions and break G4's BeamOn() into 3 pieces
    void BeamOnBeginRun( unsigned int runNumber, const char* macroFile=0, G4int n_select=-1 );
    void BeamOnDoOneArtEvent( int eventNumber );
    void BeamOnEndRun();

    void DoVisualizationFromMacro();

    std::unique_ptr<G4RunManager>         _runManager;
    std::unique_ptr<Mu2eG4ScoringManager> _scorer;

    // Do we issue warnings about multiple runs?
    bool _warnEveryNewRun;

    // Do we want to export the G4 particle data table.
    bool  _exportPDTStart;
    bool  _exportPDTEnd;

    // to be able to make StorePhysicsTable call after the event loop started
    G4VUserPhysicsList* physicsList_;
    std::string storePhysicsTablesDir_;

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
    Mu2eG4IOConfigHelper ioconf_;

    Mu2eG4PerThreadStorage perThreadStore;

    int numExcludedEvents = 0;

  }; // end G4 header

  Mu2eG4::Mu2eG4(const Parameters& pars):
    EDProducer{pars},
    conf_(pars()),
    mu2elimits_(pars().ResourceLimits()),
    trajectoryControl_(pars().TrajectoryControl()),
    multiStagePars_(pars().inputs()),
    simStage_(-1u),
    _runManager(std::make_unique<G4RunManager>()),
    _scorer(std::make_unique<Mu2eG4ScoringManager>(G4ScoringManager::GetScoringManager(),
                                                   conf_.scoring(),conf_.physics(),conf_.debug())),
    _warnEveryNewRun(pars().debug().warnEveryNewRun()),
    _exportPDTStart(pars().debug().exportPDTStart()),
    _exportPDTEnd(pars().debug().exportPDTEnd()),

    storePhysicsTablesDir_(pars().debug().storePhysicsTablesDir()),

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
    _physVolHelper(),
    _tvdOutputName(StepInstanceName::timeVD),
    timeVD_enabled_(pars().SDConfig().TimeVD().enabled()),
    _timer(std::make_unique<G4Timer>()),
    _realElapsed(0.),
    _systemElapsed(0.),
    _userElapsed(0.),
    _standardMu2eDetector((art::ServiceHandle<GeometryService>())->isStandardMu2eDetector()),
    _sensitiveDetectorHelper(pars().SDConfig()),
    ioconf_(pars(), producesCollector(), consumesCollector()),
    perThreadStore(ioconf_)
    {
      // produces() and consumes()  calls are handled by Mu2eG4IOConfigHelper

      // The string "G4Engine" is magic; see the docs for RandomNumberGenerator.
      createEngine( art::ServiceHandle<SeedService>()->getSeed(), "G4Engine");
    }

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
      validGeometryOrThrow( _rmvlevel );
    }

  }


  void Mu2eG4::initializeG4( GeometryService& geom, art::Run const& run ){

    if ( _rmvlevel > 0 ) {
      mf::LogInfo logInfo("GEOM");
      logInfo << "Initializing Geant4 for " << run.id()
              << " with verbosity " << _rmvlevel << endl;
    }


    // Create user actions and register them with G4.
    G4VUserDetectorConstruction* allMu2e = nullptr;

    // as mentioned above, we give the last element to the Master thread to setup the InstanceMap in the ctor of the SDH class

    if (_standardMu2eDetector) {
      geom.addWorldG4(*GeomHandle<Mu2eHall>());

      allMu2e =
        (new WorldMaker<Mu2eWorld>(std::make_unique<Mu2eWorld>(conf_, &(_sensitiveDetectorHelper)  ),
                                   std::make_unique<ConstructMaterials>(conf_)) );

      _originInWorld = (GeomHandle<WorldG4>())->mu2eOriginInWorld();
    }
    else {
      allMu2e =
        (new WorldMaker<Mu2eStudyWorld>(std::make_unique<Mu2eStudyWorld>(conf_, &(_sensitiveDetectorHelper) ),
                                        std::make_unique<ConstructMaterials>(conf_)) );

      // non-Mu2e detector: the system origin os set to (0.,0.,0.); do not use geometry service for that
      _originInWorld = G4ThreeVector(0.0,0.0,0.0);
    }

    preG4InitializeTasks(conf_);

    _runManager->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(allMu2e);

    physicsList_ = physicsListDecider(conf_.physics(), conf_.debug(), mu2elimits_);
    physicsList_->SetVerboseLevel(_rmvlevel);

#if G4VERSION>4106
    G4HadronicParameters::Instance()->SetVerboseLevel(_rmvlevel);
#else
    G4ParticleHPManager::GetInstance()->SetVerboseLevel(_rmvlevel);
    G4HadronicProcessStore::Instance()->SetVerbose(_rmvlevel);
#endif
    _runManager->SetUserInitialization(physicsList_);

    _scorer->initialize();

    //this is where the UserActions are instantiated
    Mu2eG4ActionInitialization* actioninit = new Mu2eG4ActionInitialization(conf_,
                                                                &_sensitiveDetectorHelper,
                                                                &perThreadStore,
                                                                &_physVolHelper,
                                                                _originInWorld
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
                                       _physVolHelper.persistentSingleStageInfo(),
                                       "");
    }
  }


  void Mu2eG4::endSubRun(art::SubRun& sr)
  {
    if(multiStagePars_.simStageOverride()) {
      const unsigned pvstage =
        writePhysicalVolumes(sr,
                             multiStagePars_.inputPhysVolumeMultiInfo(),
                             _physVolHelper.persistentSingleStageInfo(),
                             "");

      if(pvstage != simStage_) {
        throw cet::exception("BADINPUT")
          << "Mu2eG4::endSubRun() Error: inconsistent simStage: "
          <<simStage_<<" vs "<<pvstage<<"\n";
      }
    }

   _scorer->dumpInDataProduct(sr);
   _scorer->reset();

  }


  // Create one G4 event and copy its output to the art::event.
  void Mu2eG4::produce(art::Event& event) {

    perThreadStore.initializeEventInfo(&event, simStage_);

    if(multiStagePars_.updateEventLevelVolumeInfos()) {
      const unsigned pvstage =
        writePhysicalVolumes(event,
                             multiStagePars_.updateEventLevelVolumeInfos()->input,
                             _physVolHelper.persistentSingleStageInfo(),
                             multiStagePars_.updateEventLevelVolumeInfos()->outInstance);

      if(pvstage != simStage_) {
        throw cet::exception("BADINPUT")
          << "Mu2eG4::produce() Error: inconsistent simStage: "
          <<simStage_<<" vs "<<pvstage<<"\n";
      }
    }

    // Run G4 for this event and access the completed event.
    BeamOnDoOneArtEvent( event.id().event() );

    if(!perThreadStore.eventPassed()) {
      perThreadStore.clearData();
      numExcludedEvents++;
    }
    else {
      perThreadStore.putDataIntoEvent();
    }

    _runManager->TerminateOneEvent();
  }//end Mu2eG4::produce


  // Tell G4 that this run is over.
  void Mu2eG4::endRun(art::Run & run){

    BeamOnEndRun();
    if ( _rmvlevel > 0 ) {
      G4cout << "at endRun: numExcludedEvents = " << numExcludedEvents << G4endl;
    }
  }


  void Mu2eG4::endJob(){

    if ( _exportPDTEnd ) exportG4PDT( "End:" );
    _physVolHelper.endRun();
    _runManager.reset();
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

    if ( _rmvlevel > 0 ) {
      G4cout << "  Event processing inside ProcessOneEvent time summary" << G4endl;
      G4cout << "  User="  << _userElapsed
             << "s Real="  << _realElapsed
             << "s Sys="   << _systemElapsed
             << "s" << G4endl;
    }
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

DEFINE_ART_MODULE(mu2e::Mu2eG4)
