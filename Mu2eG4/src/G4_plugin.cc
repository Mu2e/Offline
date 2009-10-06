//
// A Producer Module that runs Geant4 and adds its output to the event.
// Still under development.
//
// $Id: G4_plugin.cc,v 1.2 2009/10/06 23:19:59 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/10/06 23:19:59 $
//
// Original author Rob Kutschke
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
#include <boost/shared_ptr.hpp>

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
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "Mu2eG4/inc/DetectorConstruction.hh"
#include "Mu2eG4/inc/PhysicsList.hh"
#include "Mu2eG4/inc/PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/EventAction.hh"
#include "Mu2eG4/inc/SteppingAction.hh"
#include "Mu2eG4/inc/SteppingVerbose.hh"

// ROOT includes
#include "TNtuple.h"

// This is just a placeholder for now: we need to produce something.
#include "ToyDP/inc/ToyHitCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class G4 : public edm::EDProducer {

  public:
    explicit G4(edm::ParameterSet const& pSet):
      _runManager(0),
      _genAction(0),
      _session(0),
      _visManager(0),
      _UI(0){

      // A place holder.
      produces<StepPointMCCollection>();
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
    
    G4UIsession  *_session;
    G4VisManager *_visManager;
    G4UImanager  *_UI;

    Hep3Vector _mu2eDetectorOrigin;

    // 
    string _visMacro;
    

  };
  
  // Create an instance of the run manager.
  void G4::beginJob(edm::EventSetup const&){
    _runManager = auto_ptr<Mu2eG4RunManager>(new Mu2eG4RunManager);

    edm::Service<edm::TFileService> tfs;
    
  }

  // Initialze G4.
  void G4::beginRun( edm::Run &run, edm::EventSetup const& eSetup ){

    edm::Service<GeometryService> geom;
    SimpleConfig const& config = geom->config();

    _visMacro = config.getString("graphics.macrofile","");

    static int ncalls(0);
    
    if ( ++ncalls > 1 ){
      edm::LogWarning("GEOM") 
	<< "This version of the code does not update the G4 geometry on run boundaries.";
      return;
    }

    edm::LogInfo logInfo("GEOM");
    logInfo << "Initializing Geant 4 for run: " << run.id() << endl;

    WorldMaker* allMu2e    = new WorldMaker();
    _runManager->SetUserInitialization(allMu2e);

    G4VUserPhysicsList* physics = new PhysicsList;
    _runManager->SetUserInitialization(physics);

    _genAction = new PrimaryGeneratorAction;
    _runManager->SetUserAction(_genAction);

    G4UserEventAction* event_action = new EventAction;
    _runManager->SetUserAction(event_action);
    
    SteppingAction* stepping_action = new SteppingAction;
    _runManager->SetUserAction(stepping_action);
    
    _runManager->Initialize();

    // These operations must happen after the intialize.
    // Copy some information about the G4 world to people who need it.
    Mu2eWorld const* world = allMu2e->getWorld();
    _genAction->setWorld(world);
    _mu2eDetectorOrigin = world->getMu2eDetectorOrigin();

    // Setup the graphics if requested.
    if ( _visMacro.size() > 0 ) {
      
      _UI = G4UImanager::GetUIpointer();

      _visManager = new G4VisExecutive;
      _visManager->Initialize();

      G4String command("/control/execute ");
      command += _visMacro;
      _UI->ApplyCommand( command );
      
    }
    
    // Start a run
    _runManager->BeamOnBeginRun();
    
  }

  // Create one G4 event and copy its output to the edm::event.
  void G4::produce(edm::Event& evt, edm::EventSetup const&) {

    // Create an empty output collection.
    auto_ptr<StepPointMCCollection> outputHits(new StepPointMCCollection);
    
    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<ToyGenParticleCollection> handle;
    evt.getByLabel("generate",handle);
    
    // The primary generator action needs to know about the event
    // in case it needs to get input from the event.
    _genAction->setEvent(evt);
    
    // Run G4 for this event.
    _runManager->BeamOnDoOneEvent();

    // The stuff below here should go into the end of event action.
    // Need to pass the event to that class too.

    // Access output for G4 for this event.
    G4Event const* g4event = _runManager->getCurrentEvent();
    G4HCofThisEvent* hce   = g4event->GetHCofThisEvent();

    // Get the collection ID for the straw hits.
    G4SDManager* SDman   = G4SDManager::GetSDMpointer();
    G4int colId          = SDman->GetCollectionID("StepPointG4Collection");

    if ( colId >= 0 && hce != 0 ){
      StepPointG4Collection* hits = static_cast<StepPointG4Collection*>(hce->GetHC(colId));
      G4int nHits = hits->entries();

      for (G4int i=0;i<nHits;i++) {
	StepPointG4* h = (*hits)[i];
	outputHits->push_back( h->hit() );
      }

    }

    // Should also find at the history of particles created inside G4 and copy it
    // to the edm::event.


    //_UI->ApplyCommand( "/vis/ogl/printEPS" );

    // Prompt to see graphics.
    if ( _visMacro.size() > 0 ) {
      char junk;
      cout << "Enter a character to see next event: ";
      cin >> junk;
    }
    
    evt.put(outputHits);
    
    // This deletes the object pointed to by currentEvent.
    _runManager->BeamOnEndEvent();

  }

  // Tell G4 that this run is over.
  void G4::endRun(edm::Run & run, edm::EventSetup const&){
    _runManager->BeamOnEndRun();
    if ( _visManager ) delete _visManager;
  }

  void G4::endJob(){
  }

 
} // End of namespace mu2e
 
using mu2e::G4;
DEFINE_FWK_MODULE(G4);
