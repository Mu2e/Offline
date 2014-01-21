//
// A Producer Module that runs Geant4 and adds its output to the event.
// Still under development.
//
// $Id: G4_module.cc,v 1.80 2014/01/21 06:20:41 kutschke Exp $
// $Author: kutschke $
// $Date: 2014/01/21 06:20:41 $
//
// Original author Rob Kutschke
//
//
// Notes:
// 1) According to Sunanda Banerjee, the various SetUserAction methods
//    take ownership of the object that is passed to it.  So we must
//    not delete them.
//

// Mu2e includes
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eG4/inc/Mu2eG4RunManager.hh"
#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/addPointTrajectories.hh"
#include "Mu2eG4/inc/exportG4PDT.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "Mu2eG4/inc/DetectorConstruction.hh"
#include "Mu2eG4/inc/PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/EventAction.hh"
#include "Mu2eG4/inc/SteppingAction.hh"
#include "Mu2eG4/inc/SteppingVerbose.hh"
#include "Mu2eG4/inc/StackingAction.hh"
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/postG4InitializeTasks.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/ExtMonFNALPixelSD.hh"
#include "Mu2eG4/inc/MuonMinusConversionAtRest.hh"
#include "Mu2eUtilities/inc/DiagnosticsG4.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eG4/inc/generateFieldMap.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"
#include "Mu2eG4/inc/Mu2eVisCommands.hh"

// Data products that will be produced by this module.
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
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
#include "G4VisExecutive.hh"
#include "G4Run.hh"
#include "G4Timer.hh"
#include "G4VUserPhysicsList.hh"
#include "G4RunManagerKernel.hh"
#include "G4SDManager.hh"

// ROOT includes
#include "TNtuple.h"

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

  class G4 : public art::EDProducer {

  public:
    G4(fhicl::ParameterSet const& pSet);
    // Accept compiler supplied d'tor

    virtual void produce(art::Event& e);

    virtual void beginJob();
    virtual void endJob();

    virtual void beginRun(art::Run &r);
    virtual void endRun(art::Run &);

    virtual void beginSubRun(art::SubRun &sr);

  private:
    typedef std::vector<art::InputTag> InputTags;
    typedef std::vector<std::string> Strings;

    unique_ptr<Mu2eG4RunManager> _runManager;

    // Do we issue warnings about multiple runs?
    bool _warnEveryNewRun;

    // Do we want to export the G4 particle data table.
    bool  _exportPDTStart;
    bool  _exportPDTEnd;

    PrimaryGeneratorAction* _genAction;
    TrackingAction*         _trackingAction;
    SteppingAction*         _steppingAction;
    StackingAction*         _stackingAction;

    G4UIsession  *_session;
    G4UImanager  *_UI;
    std::unique_ptr<G4VisManager> _visManager;
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
    bool   _doWriteLegacyPhysVolumeInfo;

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

    // A class to make some standard histograms.
    DiagnosticsG4 _diagnostics;

    // A parameter extracted from the geometry file at beginRun and used in produce.
    int _pointTrajectoryMinSteps;

    // Do the G4 initialization that must be done only once per job, not once per run
    void initializeG4( GeometryService& geom, art::Run const& run );

  }; // end G4 header

  G4::G4(fhicl::ParameterSet const& pSet):
    _runManager(nullptr),
    _warnEveryNewRun(pSet.get<bool>("warnEveryNewRun",false)),
    _exportPDTStart(pSet.get<bool>("exportPDTStart",false)),
    _exportPDTEnd(pSet.get<bool>("exportPDTEnd",false)),
    _genAction(nullptr),
    _trackingAction(nullptr),
    _steppingAction(nullptr),
    _stackingAction(nullptr),
    _session(nullptr),
    _UI(nullptr),
    _visManager(nullptr),
    _rmvlevel(pSet.get<int>("diagLevel",0)),
    _tmvlevel(pSet.get<int>("trackingVerbosityLevel",0)),
    _checkFieldMap(pSet.get<int>("checkFieldMap",0)),
    _visMacro(pSet.get<std::string>("visMacro","")),
    _visGUIMacro(pSet.get<std::string>("visGUIMacro","")),
    _g4Macro(pSet.get<std::string>("g4Macro","")),
    _generatorModuleLabel(pSet.get<std::string>("generatorModuleLabel", "")),
    _inputPhysVolumeMultiInfoLabel(pSet.get<string>("inputPhysVolumeMultiInfoLabel", "")),
    _doWriteLegacyPhysVolumeInfo(pSet.get<bool>("doWriteLegacyPhysVolumeInfo", true)),
    _physVolHelper(),
    _processInfo(),
    _printPhysicsProcessSummary(false),
    _simParticlePrinter(pSet.get<fhicl::ParameterSet>("SimParticlePrinter", SimParticleCollectionPrinter::defaultPSet())),
    _sensitiveDetectorHelper(pSet.get<fhicl::ParameterSet>("SDConfig", fhicl::ParameterSet())),
    _extMonFNALPixelSD(),
    _tvdOutputName(StepInstanceName::timeVD),
    _simParticleNumberOffset(pSet.get<unsigned>("simParticleNumberOffset", 0)),
    _inputSimParticles(pSet.get<std::string>("inputSimParticles", "")),
    _diagnostics(),
    _pointTrajectoryMinSteps(-1){

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
    produces<PhysicalVolumeInfoCollection,art::InRun>();
    produces<PhysicalVolumeInfoMultiCollection,art::InSubRun>();

    // The string "G4Engine" is magic; see the docs for RandomNumberGenerator.
    createEngine( art::ServiceHandle<SeedService>()->getSeed(), "G4Engine");

  } // end G4:G4(fhicl::ParameterSet const& pSet);

  // Create an instance of the run manager.
  void G4::beginJob(){
    _runManager = unique_ptr<Mu2eG4RunManager>(new Mu2eG4RunManager);
  }

  void G4::beginRun( art::Run &run){

    static int ncalls(0);
    ++ncalls;

    art::ServiceHandle<GeometryService> geom;

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
    _runManager->BeamOnBeginRun( run.id().run() );

    // Helps with indexology related to persisting G4 volume information.
    _physVolHelper.beginRun();
    _processInfo.beginRun();

    if(_doWriteLegacyPhysVolumeInfo) {
      // Add info about the G4 volumes to the run-data.
      // The framework rules requires we make a copy and add the copy.
      const PhysicalVolumeInfoCollection& vinfo = _physVolHelper.persistentInfo();
      unique_ptr<PhysicalVolumeInfoCollection> volumes(new PhysicalVolumeInfoCollection(vinfo));
      run.put(std::move(volumes));
    }

    // Some of the user actions have beginRun methods.
    GeomHandle<WorldG4>  worldGeom;
    _trackingAction->beginRun( _physVolHelper, _processInfo, worldGeom->mu2eOriginInWorld() );
    _steppingAction->beginRun( _processInfo, worldGeom->mu2eOriginInWorld() );
    _stackingAction->beginRun( worldGeom->dirtG4Ymin(), worldGeom->dirtG4Ymax() );

    // A few more things that only need to be done only once per job,
    // not once per run, but which need to be done after the call to
    // BeamOnBeginRun.

    if ( ncalls == 1 ) {

      _steppingAction->finishConstruction();

      if( _checkFieldMap>0 ) generateFieldMap(worldGeom->mu2eOriginInWorld(),_checkFieldMap);

      if ( _exportPDTStart ) exportG4PDT( "Start:" );
    }

    // Get some run-time configuration information that is stored in the geometry file.
    SimpleConfig const& config  = geom->config();
    _printPhysicsProcessSummary = config.getBool("g4.printPhysicsProcessSummary",false);

    _pointTrajectoryMinSteps = config.getInt("g4.pointTrajectoryMinSteps",5);

  }

  void G4::initializeG4( GeometryService& geom, art::Run const& run ){

    SimpleConfig const& config = geom.config();

    geom.addWorldG4();

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

    G4VUserPhysicsList* pL = physicsListDecider(config);
    pL->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(pL);

    _genAction = new PrimaryGeneratorAction();
    _runManager->SetUserAction(_genAction);

    _steppingAction = new SteppingAction(config);
    _runManager->SetUserAction(_steppingAction);

    G4UserEventAction* event_action = new EventAction(_steppingAction);
    _runManager->SetUserAction(event_action);

    _stackingAction = new StackingAction(config);
    _runManager->SetUserAction(_stackingAction);

    _trackingAction = new TrackingAction(config,_steppingAction);
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
    postG4InitializeTasks(config);
    _sensitiveDetectorHelper.registerSensitiveDetectors();
    _extMonFNALPixelSD =
      dynamic_cast<ExtMonFNALPixelSD*>(G4SDManager::GetSDMpointer()
                                       ->FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL()));


    // Setup the graphics if requested.
    if ( !_visMacro.empty() ) {

      _visManager = std::unique_ptr<G4VisManager>(new G4VisExecutive());
      _visManager->Initialize();

      ConfigFileLookupPolicy visPath;

      G4String command("/control/execute ");
      command += visPath(_visMacro);

      _UI->ApplyCommand( command );

    }

    // Book some diagnostic histograms.
    art::ServiceHandle<art::TFileService> tfs;
    _diagnostics.book("Outputs");

  } // end G4::initializeG4


  void G4::beginSubRun(art::SubRun& sr) {
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
  void G4::produce(art::Event& event) {

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

    // Some of the user actions have begin event methods. These are not G4 standards.
    _trackingAction->beginEvent(inputSimHandle, spHelper, parentHelper, *mcTrajectories );
    _genAction->setEventData(gensHandle.isValid() ? &*gensHandle : 0, genInputHits, &parentHelper);
    _steppingAction->BeginOfEvent(*tvdHits,  spHelper);

    // Connect the newly created StepPointMCCollections to their sensitive detector objects.
    _sensitiveDetectorHelper.updateSensitiveDetectors( _processInfo, spHelper);
    if(_extMonFNALPixelSD) {
      _extMonFNALPixelSD->beforeG4Event(extMonFNALHits.get(), spHelper);
    }

    // Run G4 for this event and access the completed event.
    _runManager->BeamOnDoOneEvent( event.id().event() );
    G4Event const* g4event = _runManager->getCurrentEvent();

    // Populate the output data products.
    GeomHandle<WorldG4>  world;
    GeomHandle<Mu2eBuilding>  building;
    addPointTrajectories( g4event, *pointTrajectories, spHelper, world->mu2eOriginInWorld(), _pointTrajectoryMinSteps);

    // Run self consistency checks if enabled.
    _trackingAction->endEvent(*simParticles);

    // Fill the status object.
    G4Timer const* timer = _runManager->getG4Timer();
    float cpuTime  = timer->GetSystemElapsed()+timer->GetUserElapsed();

    int status(0);
    if (  _steppingAction->nKilledStepLimit() > 0 ) status =  1;
    if (  _trackingAction->overflowSimParticles() ) status = 10;

    unique_ptr<StatusG4> g4stat(new StatusG4( status,
                                            _trackingAction->nG4Tracks(),
                                            _trackingAction->overflowSimParticles(),
                                            _steppingAction->nKilledStepLimit(),
                                            cpuTime,
                                            timer->GetRealElapsed() )
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
	  
	  if (UIE->IsGUI()) {

	    // we add a command here and initialize it (/vis/sceneHandler has to exist prior to this)
	    Mu2eVisCommandSceneHandlerDrawEvent* drEv = new Mu2eVisCommandSceneHandlerDrawEvent();
	    _visManager->RegisterMessenger(drEv); // assumes ownership;
	    // drEv->SetVisManager(_visManager.get());  
	    // vis manager pointer is static member of the drEv base class so the above is not needed

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
          UIE->SessionStart(); 
          delete UIE;

	  //If current scene is scene-0 and if scene-handler-0 has viewer-0 we
	  //will select it if not current to deal with a case which may occur
	  //e.g. in a simultaneous use of OGL & Qt

	  // basically _UI->ApplyCommand("/vis/viewer/select viewer-0"); // to have tracks drawn

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

	} // end c == 'q'

      } // end !userinput.empty()

    }   // end !_visMacro.empty()

    // This deletes the object pointed to by currentEvent.
    _runManager->BeamOnEndEvent();

  }

  // Tell G4 that this run is over.
  void G4::endRun(art::Run & run){
    _runManager->BeamOnEndRun();
  }

  void G4::endJob(){

    if ( _exportPDTEnd ) exportG4PDT( "End:" );

    // Yes, these are named endRun, but they are really endJob actions.
    _physVolHelper.endRun();
    _trackingAction->endRun();

    if ( _printPhysicsProcessSummary ){
      _processInfo.endRun();
    }

  }

} // End of namespace mu2e

using mu2e::G4;
DEFINE_ART_MODULE(G4);
