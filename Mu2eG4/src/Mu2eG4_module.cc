// A Producer Module that runs Geant4 and adds its output to the event.
//
// Original author Rob Kutschke
//
// Notes:
// 1) According to Sunanda Banerjee, the various SetUserAction methods
//    take ownership of the object that is passed to it.  So we must
//    not delete them.

// Mu2e includes
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/IMu2eG4Cut.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/addPointTrajectories.hh"
#include "Mu2eG4/inc/exportG4PDT.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eG4/inc/DetectorConstruction.hh"
#include "Mu2eG4/inc/PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/Mu2eG4SteppingAction.hh"
#include "Mu2eG4/inc/Mu2eG4StackingAction.hh"
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/postG4InitializeTasks.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/ExtMonFNALPixelSD.hh"
#include "Mu2eUtilities/inc/DiagnosticsG4.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eG4/inc/generateFieldMap.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined  G4VIS_USE_OPENGLQT ) 
#include "Mu2eG4/inc/Mu2eVisCommands.hh"
#endif

// Data products that will be produced by this module.
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"

// From art and its tool chain.
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Utilities/InputTag.h"

// Geant4 includes
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined  G4VIS_USE_OPENGLQT ) 
#include "G4VisExecutive.hh"
#endif
#include "G4Run.hh"
#include "G4Timer.hh"
#include "G4VUserPhysicsList.hh"
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

namespace mu2e {

  class Mu2eG4 : public art::EDProducer {

  public:
    Mu2eG4(fhicl::ParameterSet const& pSet);
    // Accept compiler supplied d'tor

    virtual void produce(art::Event& e) override;

    virtual void endJob() override;

    virtual void beginRun(art::Run &r) override;
    virtual void endRun(art::Run &) override;

    virtual void beginSubRun(art::SubRun &sr) override;

  private:
    fhicl::ParameterSet pset_;

    typedef std::vector<art::InputTag> InputTags;
    typedef std::vector<std::string> Strings;

    // the remnants of Mu2eG4RunManager

    // The four functions that call new G4RunManger functions and braeak the BeamOn into 4 pieces.
    void BeamOnBeginRun( unsigned int runNumber, const char* macroFile=0, G4int n_select=-1);
    void BeamOnDoOneEvent( int eventNumber );
    void BeamOnEndEvent();
    void BeamOnEndRun();

    unique_ptr<G4RunManager> _runManager;

    // Do we issue warnings about multiple runs?
    bool _warnEveryNewRun;

    // Do we want to export the G4 particle data table.
    bool  _exportPDTStart;
    bool  _exportPDTEnd;

    PrimaryGeneratorAction* _genAction;
    TrackingAction*         _trackingAction;
    Mu2eG4SteppingAction*   _steppingAction;
    Mu2eG4StackingAction*   _stackingAction;

    std::unique_ptr<IMu2eG4Cut> stackingCuts_;
    std::unique_ptr<IMu2eG4Cut> steppingCuts_;
    std::unique_ptr<IMu2eG4Cut> commonCuts_;

    G4UIsession  *_session;
    G4UImanager  *_UI;
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined  G4VIS_USE_OPENGLQT ) 
    std::unique_ptr<G4VisManager> _visManager;
#endif
    int _rmvlevel;
    int _tmvlevel;
    int _checkFieldMap;

    // Names of macro files for visualization.
    string _visMacro;  // init
    string _visGUIMacro; // end of Event GUI

    // Name of a macro file to be used for controling G4 parameters after
    // the initialization phase.
    string _g4Macro;

    art::InputTag _generatorModuleLabel;
    InputTags _genInputHitLabels;

    string _inputPhysVolumeMultiInfoLabel;

    // Helps with indexology related to persisting info about G4 volumes.
    PhysicalVolumeHelper _physVolHelper;

    // Helps with recording information about physics processes.
    PhysicsProcessInfo _processInfo;
    bool _printPhysicsProcessSummary;

    SimParticleCollectionPrinter _simParticlePrinter;

    SensitiveDetectorHelper _sensitiveDetectorHelper;
    ExtMonFNALPixelSD       *_extMonFNALPixelSD;

    // Instance name of the timeVD StepPointMC data product.
    const StepInstanceName _tvdOutputName;

    unsigned _simParticleNumberOffset;
    art::InputTag _inputSimParticles;
    art::InputTag _inputMCTrajectories;

    // A class to make some standard histograms.
    DiagnosticsG4 _diagnostics;

    // A parameter extracted from the geometry file at beginRun and used in produce.
    int _pointTrajectoryMinSteps;

    // Do the G4 initialization that must be done only once per job, not once per run
    void initializeG4( GeometryService& geom, art::Run const& run );

    std::unique_ptr<G4Timer> _timer; // local Mu2e per Geant4 event timer
    // Counters for cumulative time spent processing events by Geant4
    G4double _realElapsed;
    G4double _systemElapsed;
    G4double _userElapsed;

    // throws if obsolete config parameters are detected
    static void checkConfigRelics(const SimpleConfig& config);

  }; // end G4 header

  Mu2eG4::Mu2eG4(fhicl::ParameterSet const& pSet):
    pset_(pSet),
    _runManager(std::make_unique<G4RunManager>()),
    _warnEveryNewRun(pSet.get<bool>("debug.warnEveryNewRun",false)),
    _exportPDTStart(pSet.get<bool>("debug.exportPDTStart",false)),
    _exportPDTEnd(pSet.get<bool>("debug.exportPDTEnd",false)),
    _genAction(nullptr),
    _trackingAction(nullptr),
    _steppingAction(nullptr),
    _stackingAction(nullptr),

    stackingCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4StackingOnlyCut", fhicl::ParameterSet()))),
    steppingCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4SteppingOnlyCut", fhicl::ParameterSet()))),
    commonCuts_(createMu2eG4Cuts(pSet.get<fhicl::ParameterSet>("Mu2eG4CommonCut", fhicl::ParameterSet()))),

    _session(nullptr),
    _UI(nullptr),
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined  G4VIS_USE_OPENGLQT ) 
    _visManager(nullptr),
#endif
    // FIXME:  naming of pset parameters
    _rmvlevel(pSet.get<int>("debug.diagLevel",0)),
    _tmvlevel(pSet.get<int>("debug.trackingVerbosityLevel",0)),
    _checkFieldMap(pSet.get<int>("debug.checkFieldMap",0)),
    _visMacro(pSet.get<std::string>("visualization.macro","")),
    _visGUIMacro(pSet.get<std::string>("visualization.GUIMacro","")),
    _g4Macro(pSet.get<std::string>("g4Macro","")),
    _generatorModuleLabel(pSet.get<std::string>("generatorModuleLabel", "")),
    _inputPhysVolumeMultiInfoLabel(pSet.get<string>("inputPhysVolumeMultiInfoLabel", "")),
    _physVolHelper(),
    _processInfo(),
    _printPhysicsProcessSummary(pSet.get<bool>("printPhysicsProcessSummary",false)),
    _simParticlePrinter(pSet.get<fhicl::ParameterSet>("SimParticlePrinter", SimParticleCollectionPrinter::defaultPSet())),
    _sensitiveDetectorHelper(pSet.get<fhicl::ParameterSet>("SDConfig", fhicl::ParameterSet())),
    _extMonFNALPixelSD(),
    _tvdOutputName(StepInstanceName::timeVD),
    _simParticleNumberOffset(pSet.get<unsigned>("simParticleNumberOffset", 0)),
    _inputSimParticles(pSet.get<std::string>("inputSimParticles", "")),
    _inputMCTrajectories(pSet.get<std::string>("inputMCTrajectories", "")),
    _diagnostics(),
    _pointTrajectoryMinSteps(pSet.get<int>("g4.pointTrajectoryMinSteps",5)),
    _timer(std::make_unique<G4Timer>()),
    _realElapsed(0.),
    _systemElapsed(0.),
    _userElapsed(0.)
{

    Strings genHitsStr(pSet.get<Strings>("genInputHits", Strings()));
    for(const auto& s : genHitsStr) {
      _genInputHitLabels.emplace_back(s);
    }

    if((_generatorModuleLabel == art::InputTag()) && _genInputHitLabels.empty()) {
      throw cet::exception("CONFIG")
        << "Error: both generatorModuleLabel and genInputHits are empty - nothing to do!\n";
    }

    produces<StatusG4>();
    produces<SimParticleCollection>();

    // The main group of StepPointMCCollections.
    vector<string> const& instanceNames = _sensitiveDetectorHelper.stepInstanceNamesToBeProduced();
    for ( vector<string>::const_iterator i=instanceNames.begin();
          i != instanceNames.end(); ++i){
      produces<StepPointMCCollection>(*i);
    }

    // The timevd collection is special.
    produces<StepPointMCCollection>(_tvdOutputName.name());

    produces<PointTrajectoryCollection>();
    produces<MCTrajectoryCollection>();
    produces<ExtMonFNALSimHitCollection>();

    stackingCuts_->declareProducts(this);
    steppingCuts_->declareProducts(this);
    commonCuts_->declareProducts(this);

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

    // Helps with indexology related to persisting G4 volume information.
    _physVolHelper.beginRun();
    _processInfo.beginRun();

    // Some of the user actions have beginRun methods.
    GeomHandle<WorldG4>  worldGeom;
    _trackingAction->beginRun( _physVolHelper, _processInfo, worldGeom->mu2eOriginInWorld() );
    _steppingAction->beginRun( _processInfo, worldGeom->mu2eOriginInWorld() );

    // A few more things that only need to be done only once per job,
    // not once per run, but which need to be done after the call to
    // BeamOnBeginRun.


    if ( ncalls == 1 ) {

      _steppingAction->finishConstruction();
      stackingCuts_->finishConstruction(worldGeom->mu2eOriginInWorld());
      steppingCuts_->finishConstruction(worldGeom->mu2eOriginInWorld());
      commonCuts_->finishConstruction(worldGeom->mu2eOriginInWorld());

      if( _checkFieldMap>0 ) generateFieldMap(worldGeom->mu2eOriginInWorld(),_checkFieldMap);

      if ( _exportPDTStart ) exportG4PDT( "Start:" );
    }
  }

  void Mu2eG4::initializeG4( GeometryService& geom, art::Run const& run ){

    SimpleConfig const& config = geom.config();

    geom.addWorldG4(*GeomHandle<Mu2eHall>());

    if ( _rmvlevel > 0 ) {
      mf::LogInfo logInfo("GEOM");
      logInfo << "Initializing Geant 4 for " << run.id()
              << " with verbosity " << _rmvlevel << endl;
      logInfo << "Configured simParticleNumberOffset = "<< _simParticleNumberOffset << endl;
    }

    // Create user actions and register them with G4.

    WorldMaker<Mu2eWorld>* allMu2e    = new WorldMaker<Mu2eWorld>(std::unique_ptr<Mu2eWorld>(new Mu2eWorld(&_sensitiveDetectorHelper)));

    _runManager->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(allMu2e);

    G4VUserPhysicsList* pL = physicsListDecider(pset_);
    pL->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(pL);

    _genAction = new PrimaryGeneratorAction();
    _runManager->SetUserAction(_genAction);

    _steppingAction = new Mu2eG4SteppingAction(pset_, *steppingCuts_, *commonCuts_);
    _runManager->SetUserAction(_steppingAction);

    _stackingAction = new Mu2eG4StackingAction(pset_, *stackingCuts_, *commonCuts_);
    _runManager->SetUserAction(_stackingAction);

    _trackingAction = new TrackingAction(pset_, _steppingAction);
    _runManager->SetUserAction(_trackingAction);

    // setting tracking/stepping verbosity level; tracking manager
    // sets stepping verbosity level as well; 

    G4RunManagerKernel const * rmk = G4RunManagerKernel::GetRunManagerKernel();
    G4TrackingManager* tm  = rmk->GetTrackingManager();
    tm->SetVerboseLevel(_tmvlevel);

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

    // At this point G4 geometry and physics processes have been initialized.
    // So it is safe to modify physics processes and to compute information
    // that is derived from the G4 geometry or physics processes.

    // Mu2e specific customizations that must be done after the call to Initialize.
    postG4InitializeTasks(pset_);
    _sensitiveDetectorHelper.registerSensitiveDetectors();
    _extMonFNALPixelSD =
      dynamic_cast<ExtMonFNALPixelSD*>(G4SDManager::GetSDMpointer()
                                       ->FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL()));


#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined  G4VIS_USE_OPENGLQT ) 
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

    // Book some diagnostic histograms.
    art::ServiceHandle<art::TFileService> tfs;
    _diagnostics.book("Outputs");

  } // end G4::initializeG4


  void Mu2eG4::beginSubRun(art::SubRun& sr) {
    unique_ptr<PhysicalVolumeInfoMultiCollection> mvi(new PhysicalVolumeInfoMultiCollection());

    if(!_inputPhysVolumeMultiInfoLabel.empty()) {
      // Copy over data from the previous simulation stages
      art::Handle<PhysicalVolumeInfoMultiCollection> ih;
      sr.getByLabel(_inputPhysVolumeMultiInfoLabel, ih);
      mvi->reserve(1 + ih->size());
      mvi->insert(mvi->begin(), ih->cbegin(), ih->cend());
    }

    // Append info for the current stage
    mvi->emplace_back(std::make_pair(_simParticleNumberOffset, _physVolHelper.persistentSingleStageInfo()));

    sr.put(std::move(mvi));
  }


  // Create one G4 event and copy its output to the art::event.
  void Mu2eG4::produce(art::Event& event) {

    // Handle to the generated particles; need when building art::Ptr to a GenParticle.
    art::Handle<GenParticleCollection> gensHandle;
    if(!(_generatorModuleLabel == art::InputTag())) {
      event.getByLabel(_generatorModuleLabel, gensHandle);
    }

    // input hits from the previous simulation stage
    HitHandles genInputHits;
    for(const auto& i : _genInputHitLabels) {
      genInputHits.emplace_back(event.getValidHandle<StepPointMCCollection>(i));
    }

    art::Handle<SimParticleCollection> inputSimHandle;
    if(!(art::InputTag() == _inputSimParticles)) {
      event.getByLabel(_inputSimParticles, inputSimHandle);
      if(!inputSimHandle.isValid()) {
        throw cet::exception("CONFIG")
          << "Error retrieving inputSimParticles for "<<_inputSimParticles<<"\n";
      }
    }

    art::Handle<MCTrajectoryCollection> inputMCTracjectoryHandle;
    if(!(art::InputTag() == _inputMCTrajectories)) {
      event.getByLabel(_inputMCTrajectories, inputMCTracjectoryHandle);
      if(!inputMCTracjectoryHandle.isValid()) {
        throw cet::exception("CONFIG")
          << "Error retrieving inputMCTrajectories for "<<_inputMCTrajectories<<"\n";
      }
    }

    // ProductID for the SimParticleCollection.
    art::ProductID simPartId(getProductID<SimParticleCollection>(event));
    SimParticleHelper spHelper(_simParticleNumberOffset, simPartId, event);
    SimParticlePrimaryHelper parentHelper(event, simPartId, gensHandle);

    // Create empty data products.
    unique_ptr<SimParticleCollection>      simParticles(      new SimParticleCollection);
    unique_ptr<StepPointMCCollection>      tvdHits(           new StepPointMCCollection);
    unique_ptr<PointTrajectoryCollection>  pointTrajectories( new PointTrajectoryCollection);
    unique_ptr<MCTrajectoryCollection>     mcTrajectories(    new MCTrajectoryCollection);
    unique_ptr<ExtMonFNALSimHitCollection> extMonFNALHits(    new ExtMonFNALSimHitCollection);
    _sensitiveDetectorHelper.createProducts(event, spHelper);

    stackingCuts_->beginEvent(event, spHelper);
    steppingCuts_->beginEvent(event, spHelper);
    commonCuts_->beginEvent(event, spHelper);

    // Some of the user actions have begin event methods. These are not G4 standards.
    _trackingAction->beginEvent(inputSimHandle, inputMCTracjectoryHandle,
                                spHelper, parentHelper, *mcTrajectories );

    _genAction->setEventData(gensHandle.isValid() ? &*gensHandle : 0, genInputHits, &parentHelper);
    _steppingAction->BeginOfEvent(*tvdHits,  spHelper);

    // Connect the newly created StepPointMCCollections to their sensitive detector objects.
    _sensitiveDetectorHelper.updateSensitiveDetectors( _processInfo, spHelper);
    if(_extMonFNALPixelSD) {
      _extMonFNALPixelSD->beforeG4Event(extMonFNALHits.get(), spHelper);
    }

    // Run G4 for this event and access the completed event.
    BeamOnDoOneEvent( event.id().event() );
    G4Event const* g4event = _runManager->GetCurrentEvent();

    // Populate the output data products.
    GeomHandle<WorldG4>  world;
    addPointTrajectories( g4event, *pointTrajectories, spHelper, world->mu2eOriginInWorld(), _pointTrajectoryMinSteps);

    // Run self consistency checks if enabled.
    _trackingAction->endEvent(*simParticles);

    // Fill the status object.
    float cpuTime  = _timer->GetSystemElapsed()+_timer->GetUserElapsed();

    int status(0);
    if (  _steppingAction->nKilledStepLimit() > 0 ) status =  1;
    if (  _trackingAction->overflowSimParticles() ) status = 10;

    unique_ptr<StatusG4> g4stat(new StatusG4( status,
                                            _trackingAction->nG4Tracks(),
                                            _trackingAction->overflowSimParticles(),
                                            _steppingAction->nKilledStepLimit(),
                                            cpuTime,
                                            _timer->GetRealElapsed() )
                              );

    _diagnostics.fill( &*g4stat,
                       &*simParticles,

                       _sensitiveDetectorHelper.steps(StepInstanceName::tracker) ?
                       &_sensitiveDetectorHelper.steps(StepInstanceName::tracker).ref() : nullptr,

                       _sensitiveDetectorHelper.steps(StepInstanceName::calorimeter) ? 
                       &_sensitiveDetectorHelper.steps(StepInstanceName::calorimeter).ref() : nullptr,

                       _sensitiveDetectorHelper.steps(StepInstanceName::calorimeterRO) ?
                       &_sensitiveDetectorHelper.steps(StepInstanceName::calorimeterRO).ref() : nullptr,

                       _sensitiveDetectorHelper.steps(StepInstanceName::CRV) ?
                       &_sensitiveDetectorHelper.steps(StepInstanceName::CRV).ref() : nullptr,

                       _sensitiveDetectorHelper.steps(StepInstanceName::stoppingtarget) ?
                       &_sensitiveDetectorHelper.steps(StepInstanceName::stoppingtarget).ref() : nullptr,

                       _sensitiveDetectorHelper.steps(StepInstanceName::virtualdetector) ?
                       &_sensitiveDetectorHelper.steps(StepInstanceName::virtualdetector).ref() : nullptr,

                       _sensitiveDetectorHelper.steps(StepInstanceName::ExtMonUCITof) ?
                       &_sensitiveDetectorHelper.steps(StepInstanceName::ExtMonUCITof).ref() : nullptr,

                       &*pointTrajectories,
                       &_physVolHelper.persistentInfo() );

    _diagnostics.fillPA(_sensitiveDetectorHelper.steps(StepInstanceName::protonabsorber) ?
                        &_sensitiveDetectorHelper.steps(StepInstanceName::protonabsorber).ref() : nullptr );

    _simParticlePrinter.print(std::cout, *simParticles);

    // Add data products to the event.
    event.put(std::move(g4stat));
    event.put(std::move(simParticles));
    event.put(std::move(tvdHits),          _tvdOutputName.name()          );
    event.put(std::move(pointTrajectories));
    event.put(std::move(mcTrajectories));
    if(_extMonFNALPixelSD) {
      event.put(std::move(extMonFNALHits));
    }
    _sensitiveDetectorHelper.put(event);
    stackingCuts_->put(event);
    steppingCuts_->put(event);
    commonCuts_->put(event);

    // Pause to see graphics.
    if ( !_visMacro.empty() ){

      // Prompt to continue and wait for reply.
      cout << "Enter a character to go to the next event" << endl;
      cout << "q quits, s enters G4 interactive session, g enters a GUI session (if available)"
	   << endl;
      cout << "Once in G4 interactive session to quit it type \"exit\" "
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
	  
#if ( defined G4VIS_USE_OPENGLX || defined G4VIS_USE_OPENGL || defined  G4VIS_USE_OPENGLQT ) 

	  if (UIE->IsGUI()) {

	    // we add a command here and initialize it
	    // (/vis/sceneHandler has to exist prior to this)
	    Mu2eVisCommandSceneHandlerDrawEvent* drEv = 
              new Mu2eVisCommandSceneHandlerDrawEvent();
	    _visManager->RegisterMessenger(drEv); // assumes ownership;
	    // drEv->SetVisManager(_visManager.get());  
	    // vis manager pointer is static member of the drEv base
	    // class so the above is not needed

	    if ( !_visGUIMacro.empty() ){
	      G4String command("/control/execute ");
	      ConfigFileLookupPolicy visPath;
	      command += visPath(_visGUIMacro);
	      _UI->ApplyCommand( command );

	      cout << "In GUI interactive session use the \"Draw Current Event\" "
		   << "button in the Vis menu"
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

    }   // end !_visMacro.empty()

    // This deletes the object pointed to by currentEvent.
    BeamOnEndEvent();

  }

  // Tell G4 that this run is over.
  void Mu2eG4::endRun(art::Run & run){
    BeamOnEndRun();
  }

  void Mu2eG4::endJob(){

    if ( _exportPDTEnd ) exportG4PDT( "End:" );

    // Yes, these are named endRun, but they are really endJob actions.
    _physVolHelper.endRun();
    _trackingAction->endRun();

    if ( _printPhysicsProcessSummary ){
      _processInfo.endRun();
    }

  }

  // Do the "begin run" parts of BeamOn.
  void Mu2eG4::BeamOnBeginRun( unsigned int runNumber, const char* macroFile, G4int n_select){

    _runManager->SetRunIDCounter(runNumber);

    bool cond = _runManager->ConfirmBeamOnCondition();
    if(!cond){
      // throw here
      return;
    }

    _realElapsed   = 0.;
    _systemElapsed = 0.;
    _userElapsed   = 0.;

    //numberOfEventsToBeProcessed should be the total number of events to be processed 
    // or a large number and NOT 1 for G4 to work properly

    G4int numberOfEventsToBeProcessed = std::numeric_limits<int>::max(); // largest int for now

    _runManager->SetNumberOfEventsToBeProcessed(numberOfEventsToBeProcessed);// this would have been set by BeamOn
    _runManager->ConstructScoringWorlds();
    _runManager->RunInitialization();

    _runManager->InitializeEventLoop(numberOfEventsToBeProcessed,macroFile,n_select);

  }

  // Do the "per event" part of DoEventLoop.
  void Mu2eG4::BeamOnDoOneEvent( int eventNumber){

    // local Mu2e "per ProcessOneEvent" timer
    _timer->Start();
    _runManager->ProcessOneEvent(eventNumber);
    _timer->Stop();

    // Accumulate time spent in G4 for all events in this run.
    _realElapsed   += _timer->GetRealElapsed();
    _systemElapsed += _timer->GetSystemElapsed();
    _userElapsed   += _timer->GetUserElapsed();

  }

  void Mu2eG4::BeamOnEndEvent(){
    _runManager->TerminateOneEvent();
  }

  // Do the "end of run" parts of DoEventLoop and BeamOn.
  void Mu2eG4::BeamOnEndRun(){

    _runManager->TerminateEventLoop();

    // From G4RunManager::BeamOn.
    _runManager->RunTermination();

    G4cout << "  Event processing inside ProcessOneEvent time summary" << G4endl;
    G4cout << "  User="  << _userElapsed
           << "s Real="  << _realElapsed
           << "s Sys="   << _systemElapsed
           << "s" << G4endl;
  }

  //================================================================
  void Mu2eG4::checkConfigRelics(const SimpleConfig& config) {
    static const std::vector<std::string> keys = {

      // G4_module
      "g4.printPhysicsProcessSummary",
      "g4.pointTrajectoryMinSteps",

      // postG4InitializeTasks() call tree
      "g4.PiENuPolicy",
      "g4.PiENuPolicyVerbosity",
      "g4.minRangeCut",
      "g4.noDecay",
      "g4.doMuMinusConversionAtRest",
      "g4.useNewMuMinusAtomicCapture",

      // physicsListDecider() call tree
      "g4.physicsListName",
      "g4.useNewMuMinusAtomicCapture",
      "g4.decayMuonsWithSpin",

      // old StackingAction
      "g4.doCosmicKiller",
      "g4.cosmicKillLevel",
      "g4.cosmicVerbose",
      "g4.cosmicPcut",
      "g4.yaboveDirtYmin",
      "g4.stackPrimaryOnly",
      "g4.killLowEKine",
      "g4.killPitchToLowToStore",
      "g4.minPitch",
      "g4.stackingActionDropPDG",
      "g4.stackingActionKeepPDG",
      "g4.eKineMin",
      "g4.killLowEKinePDG",
      "g4.eKineMinPDG",

      // old TrackingAction
      "g4.particlesSizeLimit",
      "g4.mcTrajectoryMomentumCut",
      "g4.saveTrajectoryMomentumCut",
      "g4.mcTrajectoryMinSteps",
      "g4.printTrackTiming",
      "g4.trackingActionEventList",

      // old SteppingAction
      "g4.steppingActionStepsSizeLimit",
      "g4.killLowEKine",
      "g4SteppingAction.killInTheseVolumes",
      "g4SteppingAction.killerVerbose",
      "g4SteppingAction.killInHallAir",
      "g4.eKineMin",
      "g4.killLowEKinePDG",
      "g4.eKineMinPDG",
      "g4.steppingActionEventList",
      "g4.steppingActionTrackList",
      "g4.steppingActionMaxSteps",
      "g4.steppingActionMaxGlobalTime",
      "g4.steppingActionTimeVD",
      "g4.mcTrajectoryVolumes",
      "g4.mcTrajectoryVolumePtDistances",
      "g4.mcTrajectoryDefaultMinPointDistance"

    };

    std::string present;
    for(const auto k: keys) {
      if(config.hasName(k)) {
        present += k+" ";
      }
    }
    if(!present.empty()) {
      throw cet::exception("CONFIG")<<"Please use fcl to configure Mu2eG4_module. "
                                    <<"Detected obsolete SimpleConfig parameters: "<<present;
    }
  }

  //================================================================

} // End of namespace mu2e

using mu2e::Mu2eG4;
DEFINE_ART_MODULE(Mu2eG4);
