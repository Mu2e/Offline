//
// A Producer Module that runs Geant4 and adds its output to the event.
// Still under development.
//
// $Id: G4_module.cc,v 1.51 2012/06/08 22:32:18 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/06/08 22:32:18 $
//
// Original author Rob Kutschke
//
//
// Notes:
// 1) According to Sunanda Banerjee, the various SetUserAction methods
//    take ownership of the object that is passed to it.  So we must
//    not delete them.
//
// 2) Consider having a data product that is a collection type.  At
//    present there is no properly supported way for one element of the
//    collection to have an art::Ptr that refers to a different element
//    within the same collection. The issue is that we need a handle to
//    the collection before we can make the Ptrs; but we do not have
//    a handle until after the collection has become readonly.
//    The interim solution is to cast away constness.  The longer term
//    solution is to modify the post.insert() feature of data products
//    so that it can be used for a cet::map_vector<T>.
//

// Mu2e includes
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eG4/inc/Mu2eG4RunManager.hh"
#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/addPointTrajectories.hh"
#include "Mu2eG4/inc/exportG4PDT.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
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
#include "Mu2eG4/inc/MuonMinusConversionAtRest.hh"
#include "Analyses/inc/DiagnosticsG4.hh"
#include "Mu2eUtilities/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eG4/inc/generateFieldMap.hh"
#include "SeedService/inc/SeedService.hh"

// Data products that will be produced by this module.
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
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
#include "G4NistManager.hh"
#include "G4VisExecutive.hh"
#include "G4SDManager.hh"
#include "G4ParticleTable.hh"
#include "G4Run.hh"
#include "G4Timer.hh"
#include "G4VUserPhysicsList.hh"

// ROOT includes
#include "TNtuple.h"

// C++ includes.
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <iomanip>

using namespace std;
using CLHEP::Hep3Vector;

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

  private:
    auto_ptr<Mu2eG4RunManager> _runManager;

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
    G4VisManager *_visManager;
    G4UImanager  *_UI;
    int _rmvlevel;
    int _checkFieldMap;

    // Name of a macro file for visualization.
    string _visMacro;

    string _generatorModuleLabel;

    // Helps with indexology related to persisting info about G4 volumes.
    PhysicalVolumeHelper _physVolHelper;

    // Helps with recording information about physics processes.
    PhysicsProcessInfo _processInfo;
    bool _printPhysicsProcessSummary;

    // Names of output collections
    //const std::string _trackerOutputName;
    const StepInstanceName _trackerOutputName;
    const StepInstanceName      _vdOutputName;
    const StepInstanceName     _tvdOutputName;
    const StepInstanceName      _stOutputName;
    const StepInstanceName      _sbOutputName;
    const StepInstanceName    _caloOutputName;
    const StepInstanceName  _caloROOutputName;
    const StepInstanceName  _extMonFNALOutputName;
    const StepInstanceName  _extMonUCITofOutputName;
    const StepInstanceName  _ttrackerDeviceSupportOutputName;
    const StepInstanceName      _paOutputName;

    DiagnosticsG4 _diagnostics;

    // Do the G4 initialization that must be done only once per job, not once per run
    void initializeG4( GeometryService& geom, art::Run const& run );

  }; // end G4 header

  G4::G4(fhicl::ParameterSet const& pSet):
    _runManager(0),
    _warnEveryNewRun(pSet.get<bool>("warnEveryNewRun",false)),
    _exportPDTStart(pSet.get<bool>("exportPDTStart",false)),
    _exportPDTEnd(pSet.get<bool>("exportPDTEnd",false)),
    _genAction(0),
    _trackingAction(0),
    _steppingAction(0),
    _stackingAction(0),
    _session(0),
    _visManager(0),
    _UI(0),
    _rmvlevel(pSet.get<int>("diagLevel",0)),
    _checkFieldMap(pSet.get<int>("checkFieldMap",0)),
    _visMacro(pSet.get<std::string>("visMacro","")),
    _generatorModuleLabel(pSet.get<std::string>("generatorModuleLabel")),
    _physVolHelper(),
    _processInfo(),
    _printPhysicsProcessSummary(false),
    _trackerOutputName(StepInstanceName::tracker),
    _vdOutputName(StepInstanceName::virtualdetector),
    _tvdOutputName(StepInstanceName::timeVD),
    _stOutputName(StepInstanceName::stoppingtarget),
    _sbOutputName(StepInstanceName::CRV),
    _caloOutputName(StepInstanceName::calorimeter),
    _caloROOutputName(StepInstanceName::calorimeterRO),
    _extMonFNALOutputName(StepInstanceName::ExtMonFNAL),
    _extMonUCITofOutputName(StepInstanceName::ExtMonUCITof),
    _ttrackerDeviceSupportOutputName(StepInstanceName::ttrackerDS),
    _paOutputName(StepInstanceName::protonabsorber),
    _diagnostics(){

    produces<StepPointMCCollection>(_trackerOutputName.name());
    produces<StepPointMCCollection>(_vdOutputName.name());
    produces<StepPointMCCollection>(_tvdOutputName.name());
    produces<StepPointMCCollection>(_stOutputName.name());
    produces<StepPointMCCollection>(_sbOutputName.name());
    produces<StepPointMCCollection>(_caloOutputName.name());
    produces<StepPointMCCollection>(_caloROOutputName.name());
    produces<StepPointMCCollection>(_extMonFNALOutputName.name());
    produces<StepPointMCCollection>(_extMonUCITofOutputName.name());
    produces<StepPointMCCollection>(_ttrackerDeviceSupportOutputName.name());
    produces<StepPointMCCollection>(_paOutputName.name());
    produces<SimParticleCollection>();
    produces<PhysicalVolumeInfoCollection,art::InRun>();
    produces<PointTrajectoryCollection>();

    produces<StatusG4>();

    // The string "G4Engine" is magic; see the docs for RandomNumberGenerator.
    createEngine( art::ServiceHandle<SeedService>()->getSeed(), "G4Engine");

  } // end G4:G4(fhicl::ParameterSet const& pSet);

  // Create an instance of the run manager.
  void G4::beginJob(){
    _runManager = auto_ptr<Mu2eG4RunManager>(new Mu2eG4RunManager);
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

    // Add info about the G4 volumes to the run-data.
    // The framework rules requires we make a copy and add the copy.
    const PhysicalVolumeInfoCollection& vinfo = _physVolHelper.persistentInfo();
    auto_ptr<PhysicalVolumeInfoCollection> volumes(new PhysicalVolumeInfoCollection(vinfo));
    run.put(volumes);

    // Some of the user actions have beginRun methods.
    GeomHandle<WorldG4>  worldGeom;
    _trackingAction->beginRun( _physVolHelper, _processInfo, worldGeom->mu2eOriginInWorld() );
    _steppingAction->beginRun( _processInfo, worldGeom->mu2eOriginInWorld() );
    _stackingAction->beginRun( worldGeom->dirtG4Ymin(), worldGeom->dirtG4Ymax() );

    // A few more things that only need to be done only once per job, not once per run, but which need to be
    // done after the call to BeamOnBeginRun.
    if ( ncalls == 1 ) {

      _steppingAction->finishConstruction();

      if( _checkFieldMap>0 ) generateFieldMap(worldGeom->mu2eOriginInWorld(),_checkFieldMap);

      if ( _exportPDTStart ) exportG4PDT( );
    }

  }

  void G4::initializeG4( GeometryService& geom, art::Run const& run ){

    SimpleConfig const& config = geom.config();

    geom.addWorldG4();

    mf::LogInfo logInfo("GEOM");
    logInfo << "Initializing Geant 4 for " << run.id() << " with verbosity " << _rmvlevel << endl;

    // Create user actions and register them with G4.

    WorldMaker* allMu2e    = new WorldMaker();

    _runManager->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(allMu2e);

    G4VUserPhysicsList* pL = physicsListDecider(config);
    pL->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(pL);

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

    // At this point G4 geometry and physics processes have been initialized.
    // So it is safe to modify physics processes and to compute information
    // that is derived from the G4 geometry or physics processes.

    // Mu2e specific customizations that must be done after the call to Initialize.
    postG4InitializeTasks(config);

    _UI = G4UImanager::GetUIpointer();

    // Setup the graphics if requested.
    if ( !_visMacro.empty() ) {

      _visManager = new G4VisExecutive;
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

  // Create one G4 event and copy its output to the art::event.
  void G4::produce(art::Event& event) {

    // Handle to the generated particles; need when building art::Ptr to a GenParticle.
    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel(_generatorModuleLabel, gensHandle);

    // Create empty data products.
    auto_ptr<SimParticleCollection>     simParticles(      new SimParticleCollection);
    auto_ptr<StepPointMCCollection>     outputHits(        new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     vdHits(            new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     tvdHits(           new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     stHits(            new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     sbHits(            new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     caloHits(          new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     caloROHits(        new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     extMonFNALHits(    new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     extMonUCITofHits(  new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     paHits(            new StepPointMCCollection);
    auto_ptr<StepPointMCCollection>     ttrackerDeviceSupportHits(  new StepPointMCCollection);
    auto_ptr<PointTrajectoryCollection> pointTrajectories( new PointTrajectoryCollection);

    // ProductID for SimParticleCollection
    art::ProductID simPartId(getProductID<SimParticleCollection>(event));

    // Some of the user actions have begein event methods. These are not G4 standards.
    _trackingAction->beginEvent( gensHandle, simPartId, event );
    _genAction->setEvent(event);
    _steppingAction->BeginOfEvent(*tvdHits,  simPartId, event );

    // enable Sensitive Detectors to store the Framework Data Products

    G4SDManager* SDman      = G4SDManager::GetSDMpointer();

    // Get access to the master geometry system and its run time config.
    art::ServiceHandle<GeometryService> geom;
    SimpleConfig const* _config = &(geom->config());

    // Get some run-time configuration information.
    _printPhysicsProcessSummary  = _config->getBool("g4.printPhysicsProcessSummary",false);

    if ( _config->getBool("hasITracker",false) ) {
      static_cast<Mu2eSensitiveDetector*>
        (SDman->FindSensitiveDetector(SensitiveDetectorName::TrackerGas()))->
        beforeG4Event(*outputHits, _processInfo, simPartId, event );

    }else {
      static_cast<Mu2eSensitiveDetector*>
        (SDman->FindSensitiveDetector(SensitiveDetectorName::TrackerGas()))->
        beforeG4Event(*outputHits, _processInfo, simPartId, event );

      static_cast<Mu2eSensitiveDetector*>
        (SDman->FindSensitiveDetector(SensitiveDetectorName::TTrackerDeviceSupport()))->
        beforeG4Event(*ttrackerDeviceSupportHits, _processInfo, simPartId, event );
    }

    static_cast<Mu2eSensitiveDetector*>
      (SDman->FindSensitiveDetector(SensitiveDetectorName::VirtualDetector()))->
      beforeG4Event(*vdHits, _processInfo, simPartId, event );

    static_cast<Mu2eSensitiveDetector*>
      (SDman->FindSensitiveDetector(SensitiveDetectorName::StoppingTarget()))->
      beforeG4Event(*stHits, _processInfo, simPartId, event );

    static_cast<Mu2eSensitiveDetector*>
      (SDman->FindSensitiveDetector(SensitiveDetectorName::CRSScintillatorBar()))->
      beforeG4Event(*sbHits, _processInfo, simPartId, event );

    if ( geom->hasElement<Calorimeter>() ) {
      static_cast<Mu2eSensitiveDetector*>
        (SDman->FindSensitiveDetector(SensitiveDetectorName::CaloCrystal()))->
        beforeG4Event(*caloHits, _processInfo, simPartId, event );

      static_cast<Mu2eSensitiveDetector*>
        (SDman->FindSensitiveDetector(SensitiveDetectorName::CaloReadout()))->
        beforeG4Event(*caloROHits, _processInfo, simPartId, event );
    }

    static_cast<Mu2eSensitiveDetector*>
      (SDman->FindSensitiveDetector(SensitiveDetectorName::ExtMonFNAL()))->
      beforeG4Event(*extMonFNALHits, _processInfo, simPartId, event );

    static_cast<Mu2eSensitiveDetector*>
      (SDman->FindSensitiveDetector(SensitiveDetectorName::ExtMonUCITof()))->
      beforeG4Event(*extMonUCITofHits, _processInfo, simPartId, event );

    static_cast<Mu2eSensitiveDetector*>
      (SDman->FindSensitiveDetector(SensitiveDetectorName::ProtonAbsorber()))->
      beforeG4Event(*paHits, _processInfo, simPartId, event );

    // Run G4 for this event and access the completed event.
    _runManager->BeamOnDoOneEvent( event.id().event() );
    G4Event const* g4event = _runManager->getCurrentEvent();

    // Populate the output data products.
    GeomHandle<WorldG4>  world;
    GeomHandle<Mu2eBuilding>  building;
    addPointTrajectories( g4event, *pointTrajectories, building->trackerOriginInMu2e() + world->mu2eOriginInWorld());

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

    _diagnostics.fill( *g4stat,
                       *simParticles,
                       *outputHits,
                       *caloHits,
                       *caloROHits,
                       *sbHits,
                       *stHits,
                       *vdHits,
                       *extMonUCITofHits,
                       *pointTrajectories,
                       _physVolHelper.persistentInfo() );
    _diagnostics.fillPA(*paHits);

    // Add data products to the event.
    event.put(g4stat);
    event.put(simParticles);
    event.put(outputHits,        _trackerOutputName.name()     );
    event.put(vdHits,           _vdOutputName.name()           );
    event.put(tvdHits,          _tvdOutputName.name()          );
    event.put(stHits,           _stOutputName.name()           );
    event.put(sbHits,           _sbOutputName.name()           );
    event.put(caloHits,         _caloOutputName.name()         );
    event.put(caloROHits,       _caloROOutputName.name()       );
    event.put(extMonFNALHits,   _extMonFNALOutputName.name()   );
    event.put(extMonUCITofHits, _extMonUCITofOutputName.name() );
    event.put(ttrackerDeviceSupportHits, _ttrackerDeviceSupportOutputName.name() );
    event.put(paHits,           _paOutputName.name()           );
    event.put(pointTrajectories);

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
  void G4::endRun(art::Run & run){
    _runManager->BeamOnEndRun();
  }

  void G4::endJob(){

    if ( _exportPDTEnd ) exportG4PDT();

    // Yes, these are named endRun, but they are really endJob actions.
    _physVolHelper.endRun();
    _trackingAction->endRun();

    if ( _printPhysicsProcessSummary ){
      _processInfo.endRun();
    }

    delete _visManager;
  }

} // End of namespace mu2e

using mu2e::G4;
DEFINE_ART_MODULE(G4);
