//
// A Producer Module that runs Geant4 and adds its output to the event.
// ******Meant for Geant4 Studies not for Mu2e Simulations**********
//
// $Id: Mu2eG4Study_module.cc,v 1.5 2013/03/15 21:16:33 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 21:16:33 $
//
// Original author K. Genser, based on Rob's G4_module
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
#include "Mu2eG4/inc/Mu2eStudyWorld.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/addPointTrajectories.hh"
#include "Mu2eG4/inc/exportG4PDT.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "Mu2eG4/inc/DetectorConstruction.hh"
#include "Mu2eG4/inc/PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/StudyEventAction.hh"
#include "Mu2eG4/inc/StudySteppingAction.hh"
#include "Mu2eG4/inc/SteppingVerbose.hh"
#include "Mu2eG4/inc/StudyTrackingAction.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/postG4InitializeTasks.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "Mu2eG4/inc/MuonMinusConversionAtRest.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"

// Data products that will be produced by this module.
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
//#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"

// From art and its tool chain.
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// Geant4 includes
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4Run.hh"
#include "G4Timer.hh"
#include "G4VUserPhysicsList.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"

// ROOT includes
#include "TNtuple.h"

// C++ includes.
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>
#include <iomanip>

using namespace std;

namespace mu2e {

  class Mu2eG4Study : public art::EDProducer {

  public:
    Mu2eG4Study(fhicl::ParameterSet const& pSet);
    // Accept compiler supplied d'tor

    virtual void produce(art::Event& e);

    virtual void beginJob();
    virtual void endJob();

    virtual void beginRun(art::Run &r);
    virtual void endRun(art::Run &);

  private:
    unique_ptr<Mu2eG4RunManager> _runManager;

    // Do we issue warnings about multiple runs?
    bool _warnEveryNewRun;

    // Do we want to export the G4 particle data table.
    bool  _exportPDTStart;
    bool  _exportPDTEnd;

    PrimaryGeneratorAction* _genAction;
    StudyTrackingAction*    _trackingAction;
    StudySteppingAction*    _steppingAction;

    G4UIsession  *_session;
    G4UImanager  *_UI;
    std::unique_ptr<G4VisManager> _visManager;
    int _rmvlevel;
    int _tmvlevel;
    int _checkFieldMap;

    // Name of a macro file for visualization.
    string _visMacro;

    string _generatorModuleLabel;

    // Helps with indexology related to persisting info about G4 volumes.
    PhysicalVolumeHelper _physVolHelper;

    // Helps with recording information about physics processes.
    PhysicsProcessInfo _processInfo;
    bool _printPhysicsProcessSummary;

    SensitiveDetectorHelper _sensitiveDetectorHelper;

    // Instance name of the timeVD StepPointMC data product.
    const StepInstanceName _tvdOutputName;

    // Instance name of the stepperPoints StepPointMC data product.

    const StepInstanceName _steppingPointsOutputName;

    // A class to make some standard histograms.
    //    DiagnosticsG4 _diagnostics; // fixme, we may need another one for this study module

    // Do the G4 initialization that must be done only once per job, not once per run
    void initializeG4( GeometryService& geom, art::Run const& run );

  }; // end G4 header

  Mu2eG4Study::Mu2eG4Study(fhicl::ParameterSet const& pSet):
    _runManager(nullptr),
    _warnEveryNewRun(pSet.get<bool>("warnEveryNewRun",false)),
    _exportPDTStart(pSet.get<bool>("exportPDTStart",false)),
    _exportPDTEnd(pSet.get<bool>("exportPDTEnd",false)),
    _genAction(0),
    _trackingAction(0),
    _steppingAction(0),
    _session(0),
    _UI(0),
    _visManager(nullptr),
    _rmvlevel(pSet.get<int>("diagLevel",0)),
    _tmvlevel(pSet.get<int>("trackingVerbosityLevel",0)),
    _checkFieldMap(pSet.get<int>("checkFieldMap",0)),
    _visMacro(pSet.get<std::string>("visMacro","")),
    _generatorModuleLabel(pSet.get<std::string>("generatorModuleLabel")),
    _physVolHelper(),
    _processInfo(),
    _printPhysicsProcessSummary(false),
    _sensitiveDetectorHelper(pSet),
    _tvdOutputName(StepInstanceName::timeVD),
    _steppingPointsOutputName(StepInstanceName::stepper)
  {

    //    produces<StatusG4>();
    produces<SimParticleCollection>();

    // The main group of StepPointMCCollections.
    // vector<string> const& instanceNames = _sensitiveDetectorHelper.stepInstanceNamesToBeProduced();
//     for ( vector<string>::const_iterator i=instanceNames.begin();
//           i != instanceNames.end(); ++i){
//       produces<StepPointMCCollection>(*i);
//     }

    // The timevd collection is special.
    produces<StepPointMCCollection>(_tvdOutputName.name());

    // so is the stepper one
    produces<StepPointMCCollection>(_steppingPointsOutputName.name());

    //    produces<PointTrajectoryCollection>(); // may need to revisit
    produces<PhysicalVolumeInfoCollection,art::InRun>();

    // The string "G4Engine" is magic; see the docs for RandomNumberGenerator.
    createEngine( art::ServiceHandle<SeedService>()->getSeed(), "G4Engine");

  } // end Mu2eG4Study:Mu2eG4Study(fhicl::ParameterSet const& pSet);

  // Create an instance of the run manager.
  void Mu2eG4Study::beginJob(){
    _runManager = unique_ptr<Mu2eG4RunManager>(new Mu2eG4RunManager);
  }

  void Mu2eG4Study::beginRun( art::Run &run){

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

    // Add info about the G4 volumes to the run-data.
    // The framework rules requires we make a copy and add the copy.
    const PhysicalVolumeInfoCollection& vinfo = _physVolHelper.persistentInfo();
    unique_ptr<PhysicalVolumeInfoCollection> volumes(new PhysicalVolumeInfoCollection(vinfo));
    run.put(std::move(volumes));

    //  we are working in the system with the origin set to
    //  0.,0.,0. and do not use geometry service for that

    G4ThreeVector const zeroVector(0.0,0.0,0.0);

    // Some of the user actions have beginRun methods.
    _trackingAction->beginRun( _physVolHelper, _processInfo, zeroVector);
    _steppingAction->beginRun( _processInfo, zeroVector);

    // fixme may need to revisit the above User Action Clases

    // A few more things that only need to be done only once per job,
    // not once per run, but which need to be done after the call to
    // BeamOnBeginRun.

    if ( ncalls == 1 ) {
      _steppingAction->finishConstruction();
      if ( _exportPDTStart ) exportG4PDT( "Start:" );
    }

    // Get some run-time configuration information that is stored in the geometry file.
    SimpleConfig const& config  = geom->config();
    _printPhysicsProcessSummary = config.getBool("g4.printPhysicsProcessSummary",false);

  }

  void Mu2eG4Study::initializeG4( GeometryService& geom, art::Run const& run ){

    // we use GeometryService for SimpleConfig only now
    // there is still the dependence on the geometry service itself
    // via Mu2eUniverse and VolumeInfo

    SimpleConfig const& config = geom.config();

    if ( _rmvlevel > 0 ) {
      mf::LogInfo logInfo("GEOM");
      logInfo << "Initializing Geant 4 for " << run.id()
              << " with verbosity " << _rmvlevel << endl;
    }

    // Create user actions and register them with G4.

    WorldMaker<Mu2eStudyWorld>* allMu2e    = new WorldMaker<Mu2eStudyWorld>();

    _runManager->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(allMu2e);

    G4VUserPhysicsList* pL = physicsListDecider(config);
    pL->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(pL);

    _genAction = new PrimaryGeneratorAction(_generatorModuleLabel);
    _runManager->SetUserAction(_genAction);

    _steppingAction = new StudySteppingAction(config);
    _runManager->SetUserAction(_steppingAction);

    G4UserEventAction* event_action = new StudyEventAction(_steppingAction);
    _runManager->SetUserAction(event_action);

    _trackingAction = new StudyTrackingAction(config,_steppingAction);
    _runManager->SetUserAction(_trackingAction);

    // setting tracking/stepping verbosity level; tracking manager
    // sets stepping verbosity level as well; 

    G4RunManagerKernel const * rmk = G4RunManagerKernel::GetRunManagerKernel();
    G4TrackingManager* tm  = rmk->GetTrackingManager();
    tm->SetVerboseLevel(_tmvlevel);

    // Initialize G4 for this run.
    _runManager->Initialize();

    // At this point G4 geometry and physics processes have been initialized.
    // So it is safe to modify physics processes and to compute information
    // that is derived from the G4 geometry or physics processes.

    // Mu2e specific customizations that must be done after the call to Initialize.
    postG4InitializeTasks(config);

    _UI = G4UImanager::GetUIpointer();

    // Setup the graphics if requested.
    if ( !_visMacro.empty() ) {

      _visManager = std::unique_ptr<G4VisManager>(new G4VisExecutive);
      _visManager->Initialize();

      ConfigFileLookupPolicy visPath;

      G4String command("/control/execute ");
      command += visPath(_visMacro);

      _UI->ApplyCommand( command );

    }

    // Book some diagnostic histograms.
    art::ServiceHandle<art::TFileService> tfs;
    //    _diagnostics.book("Outputs");

  } // end Mu2eG4Study::initializeG4

  // Create one G4 event and copy its output to the art::event.
  void Mu2eG4Study::produce(art::Event& event) {

    // Handle to the generated particles; need when building art::Ptr to a GenParticle.
    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel(_generatorModuleLabel, gensHandle);

    // Create empty data products.
    unique_ptr<SimParticleCollection>   simParticles(new SimParticleCollection);

    unique_ptr<StepPointMCCollection> steppingPoints(new StepPointMCCollection);

    unique_ptr<StepPointMCCollection>        tvdHits(new StepPointMCCollection);

    //    unique_ptr<PointTrajectoryCollection> pointTrajectories( new PointTrajectoryCollection);

    // ProductID for the SimParticleCollection.
    art::ProductID simPartId(getProductID<SimParticleCollection>(event));

    // Some of the user actions have begein event methods. These are not G4 standards.

    _trackingAction->beginEvent(      gensHandle, simPartId, event );

    _genAction->setEvent(event);

    _steppingAction->BeginOfEvent(*tvdHits, *steppingPoints, simPartId, event );

    // Run G4 for this event and access the completed event.
    _runManager->BeamOnDoOneEvent( event.id().event() );
    //    G4Event const* g4event = _runManager->getCurrentEvent();

    // Populate the output data products.
//     GeomHandle<WorldG4>  world;
//     GeomHandle<Mu2eBuilding>  building;
    // addPointTrajectories( g4event, *pointTrajectories, building->relicMECOOriginInMu2e() + world->mu2eOriginInWorld());

    // Run self consistency checks if enabled.
    _trackingAction->endEvent(*simParticles);

    // Fill the status object.
    // G4Timer const* timer = _runManager->getG4Timer();
    // float cpuTime  = timer->GetSystemElapsed()+timer->GetUserElapsed();

//    int status(0);
//     if (  _steppingAction->nKilledStepLimit() > 0 ) status =  1;
//     if (  _trackingAction->overflowSimParticles() ) status = 10;

//     unique_ptr<StatusG4> g4stat(new StatusG4( status,
//                                             _trackingAction->nG4Tracks(),
//                                             _trackingAction->overflowSimParticles(),
//                                             _steppingAction->nKilledStepLimit(),
//                                             cpuTime,
//                                             timer->GetRealElapsed() )
//                              );

//     _diagnostics.fill( *g4stat,
//                        *simParticles,
//                        _sensitiveDetectorHelper.steps(StepInstanceName::virtualdetector).ref(),
//                        *pointTrajectories,
//                        _physVolHelper.persistentInfo() );

    // Add data products to the event.
    // event.put(std::move(g4stat));
    event.put(std::move(simParticles));
    event.put(std::move(tvdHits), _tvdOutputName.name());
    event.put(std::move(steppingPoints), _steppingPointsOutputName.name());

    //  event.put(std::move(pointTrajectories));

    // Pause to see graphics.
    if ( !_visMacro.empty() ){

      // Prompt to continue and wait for reply.
      cout << "Enter a character to go to the next event (q quits, v enters G4 interactive session)" <<
        endl;
      cout << "(Once in G4 interactive session to quit it type exit): ";
      string userinput;
      cin >> userinput;
      G4cout << userinput << G4endl;

      // Check if user is requesting an early termination of the event loop.
      if ( !userinput.empty() ){
        // Checks only the first character; we should check first non-blank.
        char c = tolower( userinput[0] );
        if ( c == 'q' ){
          throw cet::exception("CONTROL")
            << "Early end of event loop requested inside G4, \n";
        } else if ( c == 'v' ){
          G4int argc=1;
          // Cast away const-ness; required by the G4 interface ...
          char* dummy = (char *)"dummy";
          char** argv = &dummy;
          G4UIExecutive* UIE = new G4UIExecutive(argc, argv);
          UIE->SessionStart();
          delete UIE;
        }
      } // end !userinput.empty()

    }   // end !_visMacro.empty()

    // This deletes the object pointed to by currentEvent.
    _runManager->BeamOnEndEvent();

  }

  // Tell G4 that this run is over.
  void Mu2eG4Study::endRun(art::Run & run){
    _runManager->BeamOnEndRun();
  }

  void Mu2eG4Study::endJob(){

    if ( _exportPDTEnd ) exportG4PDT( "End:" );

    // Yes, these are named endRun, but they are really endJob actions.
    _physVolHelper.endRun();
    _trackingAction->endRun();

    if ( _printPhysicsProcessSummary ){
      _processInfo.endRun();
    }

  }

} // End of namespace mu2e

using mu2e::Mu2eG4Study;
DEFINE_ART_MODULE(Mu2eG4Study);
