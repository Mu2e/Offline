//
// A Producer Module that runs Geant4 and adds its output to the event.
// Still under development.
//
// $Id: G4_module.cc,v 1.16 2011/05/24 20:03:31 wb Exp $
// $Author: wb $
// $Date: 2011/05/24 20:03:31 $
//
// Original author Rob Kutschke
//
//
// Notes:
// 1) According to Sunanda Banerjee, the various SetUserAction methods
//    take ownership of the object that is passed to it.  So we must
//    not delete them.
//

// C++ includes.
#include <iostream>
#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <iomanip>

// Framework includes
#include "art/Framework/Core/Event.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"

// Geant4 includes
#include "G4UImanager.hh"
#include "G4NistManager.hh"
#include "G4VisExecutive.hh"
#include "G4SDManager.hh"
#include "G4ParticleTable.hh"
#include "G4Run.hh"
#include "G4Timer.hh"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4RunManager.hh"
#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/addPointTrajectories.hh"
#include "GeometryService/inc/GeometryService.hh"
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
#include "Mu2eG4/inc/StrawSD.hh"
#include "Mu2eG4/inc/ITGasLayerSD.hh"
#include "Mu2eG4/inc/VirtualDetectorSD.hh"
#include "Mu2eG4/inc/StoppingTargetSD.hh"
#include "Mu2eG4/inc/CRSScintillatorBarSD.hh"
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/CaloReadoutSD.hh"
#include "Mu2eG4/inc/MuonMinusConversionAtRest.hh"
#include "Mu2eG4/inc/toggleProcesses.hh"
#include "Mu2eG4/inc/DiagnosticsG4.hh"
#include "Mu2eUtilities/inc/ConfigFileLookupPolicy.hh"

// Data products that will be produced by this module.
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"

// ROOT includes
#include "TNtuple.h"

// not sure why this needs to be here; if it is above with other
// Geant4 includes a complier error occurs...

// In file included from ./MCDataProducts/inc/SimParticle.hh:22,
//              from ./MCDataProducts/inc/SimParticleCollection.hh:16,
//              from ./Mu2eG4/inc/TrackingAction.hh:22,
//              from Mu2eG4/src/G4_plugin.cc:63:
//./Mu2eUtilities/inc/PDGCode.hh:222: error: expected identifier before numeric constant
//./Mu2eUtilities/inc/PDGCode.hh:222: error: expected `}' before numeric constant
//./Mu2eUtilities/inc/PDGCode.hh:222: error: expected unqualified-id before numeric constant
//      B0 = 511 ,
#include "G4UIExecutive.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class G4 : public art::EDProducer {

  public:
    explicit G4(fhicl::ParameterSet const& pSet):
      _runManager(0),
      _genAction(0),
      _trackingAction(0),
      _steppingAction(0),
      _stackingAction(0),
      _session(0),
      _visManager(0),
      _UI(0),
      _rmvlevel(pSet.get<int>("rmvlevel",0)),
      _visMacro(pSet.get<std::string>("visMacro","")),
      _generatorModuleLabel(pSet.get<std::string>("generatorModuleLabel")),
      _physVolHelper(),
      _processInfo(),
      _printPhysicsProcessSummary(false),
      _trackerOutputName("tracker"),
      _vdOutputName("virtualdetector"),
      _stOutputName("stoppingtarget"),
      _sbOutputName("CRV"),
      _caloOutputName("calorimeter"),
      _caloROOutputName("calorimeterRO"),
      _diagnostics(){

      produces<StepPointMCCollection>(_trackerOutputName);
      produces<StepPointMCCollection>(_vdOutputName);
      produces<StepPointMCCollection>(_stOutputName);
      produces<StepPointMCCollection>(_sbOutputName);
      produces<StepPointMCCollection>(_caloOutputName);
      produces<StepPointMCCollection>(_caloROOutputName);
      produces<SimParticleCollection>();
      produces<PhysicalVolumeInfoCollection,art::InRun>();
      produces<PointTrajectoryCollection>();

      produces<StatusG4>();

      // The string "G4Engine" is magic; see the docs for RandomNumberGenerator.
      createEngine( get_seed_value(pSet), "G4Engine");

    }

    virtual ~G4() {
      // See note 1.
    }

    virtual void produce(art::Event& e);

    virtual void beginJob();
    virtual void endJob();

    virtual void beginRun(art::Run &r);
    virtual void endRun(art::Run &);

  private:
    auto_ptr<Mu2eG4RunManager> _runManager;

    PrimaryGeneratorAction* _genAction;
    TrackingAction*         _trackingAction;
    SteppingAction*         _steppingAction;
    StackingAction*         _stackingAction;

    G4UIsession  *_session;
    G4VisManager *_visManager;
    G4UImanager  *_UI;
    int _rmvlevel;


    // Position, in G4 world coord, of (0,0,0) of the mu2e coordinate system.
    CLHEP::Hep3Vector _mu2eOrigin;

    // Position, in G4 world coord, of (0,0,0) of the detector coordinate system.
    CLHEP::Hep3Vector _mu2eDetectorOrigin;

    // Name of a macro file for visualization.
    string _visMacro;

    string _generatorModuleLabel;

    // Helps with indexology related to persisting info about G4 volumes.
    PhysicalVolumeHelper _physVolHelper;

    // Helps with recording information about physics processes.
    PhysicsProcessInfo _processInfo;
    bool _printPhysicsProcessSummary;

    // Names of output collections
    const std::string _trackerOutputName;
    const std::string      _vdOutputName;
    const std::string      _stOutputName;
    const std::string      _sbOutputName;
    const std::string    _caloOutputName;
    const std::string  _caloROOutputName;

    DiagnosticsG4 _diagnostics;

  };

  // Create an instance of the run manager.
  void G4::beginJob(){
    _runManager = auto_ptr<Mu2eG4RunManager>(new Mu2eG4RunManager);

    // If you want job scope histograms.
    art::ServiceHandle<art::TFileService> tfs;
    _diagnostics.beginJob();

  }

  // Initialze G4.
  void G4::beginRun( art::Run &run){

    art::ServiceHandle<GeometryService> geom;
    SimpleConfig const& config = geom->config();

    static int ncalls(0);

    if ( ++ncalls > 1 ){
      mf::LogWarning("GEOM")
        << "This version of the code does not update the G4 geometry on run boundaries.";
      return;
    }

    mf::LogInfo logInfo("GEOM");
    logInfo << "Initializing Geant 4 for run: " << run.id() << endl;

    // Create user actions and register them with G4.

    WorldMaker* allMu2e    = new WorldMaker();

    _runManager->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(allMu2e);

    _runManager->SetUserInitialization(physicsListDecider(config));

    _genAction = new PrimaryGeneratorAction(_generatorModuleLabel);
    _runManager->SetUserAction(_genAction);

    _steppingAction = new SteppingAction(config);
    _runManager->SetUserAction(_steppingAction);

    G4UserEventAction* event_action = new EventAction(_steppingAction);
    _runManager->SetUserAction(event_action);

    _stackingAction = new StackingAction(config);
    _runManager->SetUserAction(_stackingAction);

    _trackingAction = new TrackingAction(config,_steppingAction);
    _runManager->SetUserAction(_trackingAction);

    // Initialize G4 for this run.
    _runManager->Initialize();

    // Switch off the decay of some particles
    switchDecayOff(config);

    // add user processes
    addUserProcesses(config);

    // At this point G4 geometry has been initialized.  So it is safe to initialize
    // objects that depend on G4 geometry.

    // Copy some information about the G4 world to people who need it.
    Mu2eWorld const* world = allMu2e->getWorld();
    _mu2eOrigin            = world->getMu2eOrigin();
    _mu2eDetectorOrigin    = world->getMu2eDetectorOrigin();

    // Setup the graphics if requested.
    if ( !_visMacro.empty() ) {

      _UI = G4UImanager::GetUIpointer();

      _visManager = new G4VisExecutive;
      _visManager->Initialize();

      ConfigFileLookupPolicy visPath;

      G4String command("/control/execute ");
      command += visPath(_visMacro);

      _UI->ApplyCommand( command );

    }

    // Start a run
    _runManager->BeamOnBeginRun( run.id().run() );

    // Helps with indexology related to persisting G4 volume information.
    _physVolHelper.beginRun();
    _processInfo.beginRun();

    // Add info about the G4 volumes to the run-data.
    // The framework rules requires we make a copy and add the copy.
    const PhysicalVolumeInfoCollection& vinfo = _physVolHelper.persistentInfo();
    auto_ptr<PhysicalVolumeInfoCollection> volumes(new PhysicalVolumeInfoCollection(vinfo));
    run.put(volumes);

    // Some of the user actions have beginRun methods.
    _genAction->setWorld(world);
    _trackingAction->beginRun( _physVolHelper, _processInfo, _mu2eOrigin );
    _steppingAction->beginRun( );
    _stackingAction->beginRun( world->getDirtG4Ymin(), world->getDirtG4Ymax() );

    _diagnostics.beginRun( run, _physVolHelper );

  }

  // Create one G4 event and copy its output to the art::event.
  void G4::produce(art::Event& event) {

    // Create empty data products.
    auto_ptr<SimParticleCollection>     simParticles(      new SimParticleCollection);
    auto_ptr<StepPointMCCollection>     outputHits(        new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     vdHits(            new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     stHits(            new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     sbHits(            new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     caloHits(          new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     caloROHits(        new StepPointMCCollection);
    auto_ptr<PointTrajectoryCollection> pointTrajectories( new PointTrajectoryCollection);

     // Some of the user actions have begein event methods. These are not G4 standards.
    _trackingAction->beginEvent();
    _genAction->setEvent(event);

    // enable Sensitive Detectors to store the Framework Data Products

    G4SDManager* SDman      = G4SDManager::GetSDMpointer();

    // Get access to the master geometry system and its run time config.
    art::ServiceHandle<GeometryService> geom;
    SimpleConfig const* _config = &(geom->config());

    // Get some run-time configuration information.
    _printPhysicsProcessSummary  = _config->getBool("g4.printPhysicsProcessSummary",false);

    if ( _config->getBool("hasITracker",false) ) {
            static_cast<ITGasLayerSD*>
            (SDman->FindSensitiveDetector(SensitiveDetectorName::ItrackerGasVolume()))->
            beforeG4Event(*outputHits, _processInfo);

    }else {
            static_cast<StrawSD*>
            (SDman->FindSensitiveDetector(SensitiveDetectorName::StrawGasVolume()))->
             beforeG4Event(*outputHits, _processInfo);
    }

    static_cast<VirtualDetectorSD*>
      (SDman->FindSensitiveDetector(SensitiveDetectorName::VirtualDetector()))->
      beforeG4Event(*vdHits, _processInfo);

    static_cast<StoppingTargetSD*>
      (SDman->FindSensitiveDetector(SensitiveDetectorName::StoppingTarget()))->
      beforeG4Event(*stHits, _processInfo);

    static_cast<CRSScintillatorBarSD*>
      (SDman->FindSensitiveDetector(SensitiveDetectorName::CRSScintillatorBar()))->
      beforeG4Event(*sbHits, _processInfo);

    static_cast<CaloCrystalSD*>
      (SDman->FindSensitiveDetector(SensitiveDetectorName::CaloCrystal()))->
      beforeG4Event(*caloHits, _processInfo);

    static_cast<CaloReadoutSD*>
      (SDman->FindSensitiveDetector(SensitiveDetectorName::CaloReadout()))->
      beforeG4Event(*caloROHits, _processInfo);

    // Run G4 for this event and access the completed event.
    _runManager->BeamOnDoOneEvent( event.id().event() );
    G4Event const* g4event = _runManager->getCurrentEvent();

    // Populate the output data products.
    addPointTrajectories( g4event, *pointTrajectories, _mu2eDetectorOrigin);

    // Run self consistency checks if enabled.
    _trackingAction->endEvent(*simParticles);

    // Fill the status object.
    G4Timer const* timer = _runManager->getG4Timer();
    float cpuTime  = timer->GetSystemElapsed()+timer->GetUserElapsed();

    int status(0);
    if (  _steppingAction->nKilledStepLimit() > 0 ) status =  1;
    if (  _trackingAction->overflowSimParticles() ) status = 10;

    auto_ptr<StatusG4> g4stat(new StatusG4( status,
                                            _trackingAction->nG4Tracks(),
                                            _trackingAction->overflowSimParticles(),
                                            _steppingAction->nKilledStepLimit(),
                                            cpuTime,
                                            timer->GetRealElapsed() )
                              );

    _diagnostics.analyze( *g4stat,
                          *simParticles,
                          *outputHits,
                          *caloHits,
                          *caloROHits,
                          *sbHits,
                          *stHits,
                          *vdHits,
                          *pointTrajectories,
                          _physVolHelper);

    // Add data products to the event.
    event.put(g4stat);
    event.put(outputHits,_trackerOutputName);
    event.put(vdHits,_vdOutputName);
    event.put(stHits,_stOutputName);
    event.put(sbHits,_sbOutputName);
    event.put(simParticles);
    event.put(caloHits,_caloOutputName);
    event.put(caloROHits,_caloROOutputName);
    event.put(pointTrajectories);


    // Pause to see graphics.
    if ( !_visMacro.empty() ){

      // We need to add a command here to flush the graphics to the screen
      // The only way that I know how does a full redraw ...

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
      _UI->ApplyCommand("/vis/scene/endOfEventAction refresh");

    }   // end !_visMacro.empty()


    // This deletes the object pointed to by currentEvent.
    _runManager->BeamOnEndEvent();

  }

  // Tell G4 that this run is over.
  void G4::endRun(art::Run & run){

    _diagnostics.endRun(run);

    _runManager->BeamOnEndRun();
    _physVolHelper.endRun();
    _trackingAction->endRun();

    if ( _printPhysicsProcessSummary ){
      _processInfo.endRun();
    }

    delete _visManager;
  }

  void G4::endJob(){
    _diagnostics.endJob();
  }



} // End of namespace mu2e

using mu2e::G4;
DEFINE_ART_MODULE(G4);
