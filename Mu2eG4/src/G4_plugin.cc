//
// A Producer Module that runs Geant4 and adds its output to the event.
// Still under development.
//
// $Id: G4_plugin.cc,v 1.19 2010/05/18 21:16:14 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/05/18 21:16:14 $
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
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

// Geant4 includes
#include "G4UImanager.hh"
#include "G4NistManager.hh"
#include "G4VisExecutive.hh"
#include "G4SDManager.hh"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4RunManager.hh"
#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/StepPointG4.hh"
#include "Mu2eG4/inc/addStepPointMCs.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "Mu2eG4/inc/DetectorConstruction.hh"
#include "Mu2eG4/inc/PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/EventAction.hh"
#include "Mu2eG4/inc/SteppingAction.hh"
#include "Mu2eG4/inc/SteppingVerbose.hh"
#include "Mu2eG4/inc/StackingAction.hh"
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/physicsListDecider.hh"

#include "ITrackerGeom/inc/ITracker.hh"

// Data products that will be produced by this module.
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"

// ROOT includes
#include "TNtuple.h"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class G4 : public edm::EDProducer {

  public:
    explicit G4(edm::ParameterSet const& pSet):
      _runManager(0),
      _genAction(0),
      _trackingAction(0),
      _session(0),
      _visManager(0),
      _UI(0),
      _visMacro(pSet.getUntrackedParameter<std::string>("visMacro","")),
      _generatorModuleLabel(pSet.getParameter<std::string>("generatorModuleLabel")){

      produces<StepPointMCCollection>();
      produces<SimParticleCollection>();
      produces<PhysicalVolumeInfoCollection,edm::InRun>();

    }

    virtual ~G4() { 
      // Must not delete the pointers handed to G4.
      // G4 takes over the lifetime of these objects.
    }

    virtual void produce(edm::Event& e, edm::EventSetup const& c);
    
    virtual void beginJob(edm::EventSetup const&);
    virtual void endJob();
 
    virtual void beginRun(edm::Run &r, edm::EventSetup const& eSetup );
    virtual void endRun(edm::Run &, edm::EventSetup const&);

    static void fillDescription(edm::ParameterSetDescription& iDesc,
                                string const& moduleLabel) {
      iDesc.setAllowAnything();
    }
    
  private:
    auto_ptr<Mu2eG4RunManager> _runManager;

    PrimaryGeneratorAction* _genAction;
    TrackingAction*         _trackingAction;
    
    G4UIsession  *_session;
    G4VisManager *_visManager;
    G4UImanager  *_UI;

    // Position, in G4 world coord, of (0,0,0) of the mu2e coordinate system.
    CLHEP::Hep3Vector _mu2eOrigin;

    // Position, in G4 world coord, of (0,0,0) of the detector coordinate system.
    CLHEP::Hep3Vector _mu2eDetectorOrigin;

    // Name of a macro file for visualization.
    string _visMacro;

    string _generatorModuleLabel;

    // Helps with indexology related to persisting info about G4 volumes.
    PhysicalVolumeHelper _physVolHelper;

  };
  
  // Create an instance of the run manager.
  void G4::beginJob(edm::EventSetup const&){
    _runManager = auto_ptr<Mu2eG4RunManager>(new Mu2eG4RunManager);

    // If you want job scope histograms.
    edm::Service<edm::TFileService> tfs;
    
  }

  // Initialze G4.
  void G4::beginRun( edm::Run &run, edm::EventSetup const& eSetup ){

    edm::Service<GeometryService> geom;
    SimpleConfig const& config = geom->config();

    static int ncalls(0);
    
    if ( ++ncalls > 1 ){
      edm::LogWarning("GEOM") 
        << "This version of the code does not update the G4 geometry on run boundaries.";
      return;
    }

    edm::LogInfo logInfo("GEOM");
    logInfo << "Initializing Geant 4 for run: " << run.id() << endl;

    // Create user actions and register them with G4.

    WorldMaker* allMu2e    = new WorldMaker();
    _runManager->SetUserInitialization(allMu2e);

    _runManager->SetUserInitialization(physicsListDecider(config));

    _genAction = new PrimaryGeneratorAction(_generatorModuleLabel);
    _runManager->SetUserAction(_genAction);

    G4UserEventAction* event_action = new EventAction;
    _runManager->SetUserAction(event_action);
    
    SteppingAction* stepping_action = new SteppingAction(config);
    _runManager->SetUserAction(stepping_action);

    StackingAction* stacking_action = new StackingAction(config);
    _runManager->SetUserAction(stacking_action);

    _trackingAction = new TrackingAction(config);
    _runManager->SetUserAction(_trackingAction);

    // Initialize G4 for this run.
    _runManager->Initialize();

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

      G4String command("/control/execute ");
      command += _visMacro;
      _UI->ApplyCommand( command );
      
    }
    
    // Start a run
    _runManager->BeamOnBeginRun();

    // Helps with indexology related to persisting G4 volume information.
    _physVolHelper.beginRun();

    // Add info about the G4 volumes to the run-data.  
    // The framework rules requires we make a copy and add the copy.
    const PhysicalVolumeInfoCollection& vinfo = _physVolHelper.persistentInfo();
    auto_ptr<PhysicalVolumeInfoCollection> volumes(new PhysicalVolumeInfoCollection(vinfo));
    run.put(volumes);

    // Some of the user actions have beginRun methods.
    _genAction->setWorld(world);
    _trackingAction->beginRun( _physVolHelper, _mu2eOrigin );

  }

  // Create one G4 event and copy its output to the edm::event.
  void G4::produce(edm::Event& event, edm::EventSetup const&) {

    // Create empty data products.
    auto_ptr<StepPointMCCollection> outputHits(new StepPointMCCollection);
    auto_ptr<SimParticleCollection> simParticles(new SimParticleCollection);

    // Some of the user actions have begein event methods. These are not G4 standards.
    _trackingAction->beginEvent();
    _genAction->setEvent(event);
    
    // Run G4 for this event and access the completed event.
    _runManager->BeamOnDoOneEvent();
    G4Event const* g4event = _runManager->getCurrentEvent();

    // Populate the output data products.
    addStepPointMCs( g4event, *outputHits);
    _trackingAction->endEvent( *simParticles );

    event.put(outputHits);
    event.put(simParticles);

    // Pause to see graphics. 
    if ( _visMacro.size() > 0 ) {

      // Prompt to continue and wait for reply.
      cout << "Enter a character to see next event: "; 
      string junk;
      cin >> junk;
      
      _UI->ApplyCommand( "/vis/viewer/refresh"); 

      // Check if user is requesting an early termination of the event loop.
      if ( !junk.empty() ){

        // Checks only the first character; we should check first non-blank.
        char c = tolower( junk[0] );
        if ( c == 'q' ){
          throw cms::Exception("CONTROL")
            << "Early end of event loop requested inside G4, \n";
        }
      }
    }

    // This deletes the object pointed to by currentEvent.
    _runManager->BeamOnEndEvent();

  }

  // Tell G4 that this run is over.
  void G4::endRun(edm::Run & run, edm::EventSetup const&){
    _runManager->BeamOnEndRun();
    _physVolHelper.endRun();

    if ( _visManager ) delete _visManager;
  }

  void G4::endJob(){
  }
 
} // End of namespace mu2e

using mu2e::G4;
DEFINE_FWK_MODULE(G4);
