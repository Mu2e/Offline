//
// A Producer Module that runs Geant4 and adds its output to the event.
// Still under development.
//
// $Id: G4_plugin.cc,v 1.27 2010/09/20 02:57:05 logash Exp $
// $Author: logash $ 
// $Date: 2010/09/20 02:57:05 $
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
#include "G4ParticleTable.hh"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4RunManager.hh"
#include "Mu2eG4/inc/WorldMaker.hh"
#include "Mu2eG4/inc/Mu2eWorld.hh"
#include "Mu2eG4/inc/StepPointG4.hh"
#include "Mu2eG4/inc/addStepPointMCs.hh"
#include "Mu2eG4/inc/addVirtualDetectorPoints.hh"
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

#include "Mu2eG4/inc/VirtualDetectorSD.hh"
#include "Mu2eG4/inc/StrawSD.hh"

#include "ITrackerGeom/inc/ITracker.hh"

// Data products that will be produced by this module.
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"

// ROOT includes
#include "TNtuple.h"

// not sure why this needs to be here; if it is above with other
// Geant4 includes a complier error occurs...

// In file included from ./ToyDP/inc/SimParticle.hh:22,
//              from ./ToyDP/inc/SimParticleCollection.hh:16,
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

  class G4 : public edm::EDProducer {

  public:
    explicit G4(edm::ParameterSet const& pSet):
      _runManager(0),
      _genAction(0),
      _trackingAction(0),
      _session(0),
      _visManager(0),
      _UI(0),
      _rmvlevel(pSet.getUntrackedParameter<int>("rmvlevel",0)),
      _visMacro(pSet.getUntrackedParameter<std::string>("visMacro","")),
      _generatorModuleLabel(pSet.getParameter<std::string>("generatorModuleLabel")),
      _trackerOutputName("tracker"),
      _stepsSizeLimit(pSet.getUntrackedParameter<int>("stepsSizeLimit",10000)),
      _particlesSizeLimit(pSet.getUntrackedParameter<int>("particlesSizeLimit",10000)),
      _vdOutputName("virtualdetector") {

      produces<StepPointMCCollection>(_trackerOutputName);
      produces<StepPointMCCollection>(_vdOutputName);
      produces<SimParticleCollection>();
      produces<PhysicalVolumeInfoCollection,edm::InRun>();

      // The string "G4Engine" is magic; see the docs for RandomNumberGeneratorService.
      createEngine( get_seed_value(pSet), "G4Engine");

    }

    virtual ~G4() { 
      // See note 1.
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
    
    void switchDecayOff(const SimpleConfig&);

  private:
    auto_ptr<Mu2eG4RunManager> _runManager;

    PrimaryGeneratorAction* _genAction;
    TrackingAction*       _trackingAction;
    
    G4UIsession  *_session;
    G4VisManager *_visManager;
    G4UImanager  *_UI;
    int _rmvlevel;

    // Limits for collections size
    int _stepsSizeLimit;
    int _particlesSizeLimit;

    // Position, in G4 world coord, of (0,0,0) of the mu2e coordinate system.
    CLHEP::Hep3Vector _mu2eOrigin;

    // Position, in G4 world coord, of (0,0,0) of the detector coordinate system.
    CLHEP::Hep3Vector _mu2eDetectorOrigin;

    // Name of a macro file for visualization.
    string _visMacro;

    string _generatorModuleLabel;

    // Helps with indexology related to persisting info about G4 volumes.
    PhysicalVolumeHelper _physVolHelper;

    // Names of output collections
    const std::string _trackerOutputName;
    const std::string      _vdOutputName;

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

    _runManager->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(allMu2e);

    _runManager->SetUserInitialization(physicsListDecider(config));

    _genAction = new PrimaryGeneratorAction(_generatorModuleLabel);
    _runManager->SetUserAction(_genAction);

    SteppingAction* stepping_action = new SteppingAction(config);
    _runManager->SetUserAction(stepping_action);

    G4UserEventAction* event_action = new EventAction(stepping_action);
    _runManager->SetUserAction(event_action);
    
    StackingAction* stacking_action = new StackingAction(config);
    _runManager->SetUserAction(stacking_action);

    _trackingAction = new TrackingAction(config);
    _runManager->SetUserAction(_trackingAction);

    // Initialize G4 for this run.
    _runManager->Initialize();

    // Switch off the decay of some particles
    switchDecayOff(config);

    // At this point G4 geometry has been initialized.  So it is safe to initialize
    // objects that depend on G4 geometry.

    // Copy some information about the G4 world to people who need it.
    Mu2eWorld const* world = allMu2e->getWorld();
    _mu2eOrigin          = world->getMu2eOrigin();
    _mu2eDetectorOrigin    = world->getMu2eDetectorOrigin();

    // Limit size of output collections
    StrawSD::setSizeLimit(_stepsSizeLimit);
    VirtualDetectorSD::setSizeLimit(_stepsSizeLimit);
    TrackingAction::setSizeLimit(_particlesSizeLimit);

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
    stepping_action->beginRun();

  }

  // Create one G4 event and copy its output to the edm::event.
  void G4::produce(edm::Event& event, edm::EventSetup const&) {

    // Create empty data products.
    auto_ptr<StepPointMCCollection> outputHits(new StepPointMCCollection);
    auto_ptr<SimParticleCollection> simParticles(new SimParticleCollection);
    auto_ptr<StepPointMCCollection> vdHits(new StepPointMCCollection);

    // Some of the user actions have begein event methods. These are not G4 standards.
    _trackingAction->beginEvent();
    _genAction->setEvent(event);
    
    // Run G4 for this event and access the completed event.
    _runManager->BeamOnDoOneEvent();
    G4Event const* g4event = _runManager->getCurrentEvent();

    // Populate the output data products.
    addStepPointMCs( g4event, *outputHits);
    addVirtualDetectorPoints( g4event, *vdHits);
    _trackingAction->endEvent( *simParticles );

    event.put(outputHits,_trackerOutputName);
    event.put(vdHits,_vdOutputName);
    event.put(simParticles);
    
    //     // Pause to see graphics. 
    //     if ( _visMacro.size() > 0 ) {

    //       _UI->ApplyCommand( "/vis/scene/endOfEventAction refresh");

    //       // Prompt to continue and wait for reply.
    //       cout << "Enter a character to see next event: "; 
    //       string junk;
    //       cin >> junk;
      
    //       // Check if user is requesting an early termination of the event loop.
    //       if ( !junk.empty() ){

    //       // Checks only the first character; we should check first non-blank.
    //       char c = tolower( junk[0] );
    //       if ( c == 'q' ){
    //         throw cms::Exception("CONTROL")
    //           << "Early end of event loop requested inside G4, \n";
    //       }
    //       }
    //     }

    // Pause to see graphics. 
    if ( _visMacro.size() > 0 ) {

      _UI->ApplyCommand( "/vis/scene/endOfEventAction refresh");

      // Prompt to continue and wait for reply.
      cout << "Enter a character to go to the next event (q quits, v enters G4 interactive session)" << endl;
      cout << "(Once in G4 interactive session to quit it type exit): ";
      string userinput;
      cin >> userinput;
      G4cout << userinput << G4endl;
      
      // Check if user is requesting an early termination of the event loop.
      if ( !userinput.empty() ){
        // Checks only the first character; we should check first non-blank.
        char c = tolower( userinput[0] );
        if ( c == 'q' ){
          throw cms::Exception("CONTROL")
            << "Early end of event loop requested inside G4, \n";
        } else if ( c == 'v' ){
          G4int argc=1;
          char* dummy = "dummy";
          char** argv = &dummy;
          G4UIExecutive* UIE = new G4UIExecutive(argc, argv);
          UIE->SessionStart();
          delete UIE;
        }
      }
      //_UI->ApplyCommand(userinput); 
      _UI->ApplyCommand("/vis/scene/endOfEventAction refresh");
      //_UI->ApplyCommand( "/vis/viewer/refresh"); 
      //_UI->ApplyCommand( "/vis/viewer/flush"); 

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

  void G4::switchDecayOff(const SimpleConfig& config) {

    // Read list of particles for which the decay should be switched off
    vector<int> plist;
    if( ! config.hasName("g4.noDecay") ) return;
    config.getVectorInt("g4.noDecay",plist);
    
    G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();
    for( int i=0; i<plist.size(); ++i ) {
      int pdg = plist[i];
      G4ParticleDefinition* particle = theParticleTable->FindParticle(pdg);
      if( particle==0 ) {
	cout << "SwitchDecayOff: cannot find particle pdgId=" << pdg << endl;
      } else {
	G4ProcessManager* pmanager = particle->GetProcessManager();
	G4ProcessVector * pVector  = pmanager->GetProcessList();
	G4VProcess *decayProcess = 0;
	for( G4int j=0; j<pmanager->GetProcessListLength(); j++ ) {
	  if( (*pVector)[j]->GetProcessName() == "Decay" ) {
	    decayProcess = (*pVector)[j];
	    break;
	  }
	}
	if( decayProcess==0 ) {
	  cout << "SwitchDecayOff: cannot find decay process for particle pdgId=" << pdg 
	       << " (" << particle->GetParticleName() << ")" << endl;
	} else {
	  pmanager->RemoveProcess(decayProcess);
	  cout << "SwitchDecayOff: decay process is removed for particle pdgId=" << pdg 
	       << " (" << particle->GetParticleName() << ")" << endl;
	}
	cout << "SwitchDecayOff: list of processes defined for particle pdgId=" << pdg 
	     << " (" << particle->GetParticleName() << "):" << endl;
	for( G4int j=0; j<pmanager->GetProcessListLength(); j++ ) 
	  cout << (*pVector)[j]->GetProcessName() << endl;
      }
    }

  }
 
} // End of namespace mu2e

using mu2e::G4;
DEFINE_FWK_MODULE(G4);
