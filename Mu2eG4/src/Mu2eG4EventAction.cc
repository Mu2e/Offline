//
// G4 begin and end of event actions for Mu2e.
//
// Author: Lisa Goodeough
// Date: 2017/05/04
//

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4EventAction.hh"
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/Mu2eG4SteppingAction.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"
#include "Mu2eG4/inc/IMu2eG4Cut.hh"
#include "MCDataProducts/inc/StatusG4.hh"

//G4 includes
#include "G4Timer.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"

//art includes
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Principal/Event.h"

//C++ includes
#include <iostream>

using namespace std;

namespace {
  // Because a Mu2eG4EventAction object can be created from multiple
  // threads and because it uses configuration validation, which is
  // not thread-safe, there exists the possibility of data races.
  //
  // To avoid this, we create a global configuration struct, which
  // does all of the thread-unsafe manipulations upon construction
  // (namely adjusting FHiCL's table-member and name-stack
  // registries). We then copy the struct in the body of the c'tor,
  // where they copy operation does not use the name-stack registry.
  //
  // NB - This is an expert-only workaround which should arguably go
  //      away if fhiclcpp decides to adopt thread-local registries
  //      for configuration validation.
  mu2e::SimParticleCollectionPrinter::Config const global_c;
}

namespace mu2e {

  Mu2eG4EventAction::Mu2eG4EventAction(const Mu2eG4Config::Top& conf,
                                       TrackingAction* tracking_action,
                                       Mu2eG4SteppingAction* stepping_action,
                                       SensitiveDetectorHelper* sensitive_detectorhelper,
                                       IMu2eG4Cut& stacking_cuts,
                                       IMu2eG4Cut& stepping_cuts,
                                       IMu2eG4Cut& common_cuts,
                                       Mu2eG4PerThreadStorage* pts,
                                       PhysicsProcessInfo* phys_process_info,
                                       const CLHEP::Hep3Vector& origin_in_world)
  :
  G4UserEventAction(),
    perThreadObjects_(pts),
    simParticlePrinter_(),
    multiStagePars_(conf),
    _trackingAction(tracking_action),
    _steppingAction(stepping_action),
    _sensitiveDetectorHelper(sensitive_detectorhelper),
    _stackingCuts(&stacking_cuts),
    _steppingCuts(&stepping_cuts),
    _commonCuts(&common_cuts),
    _originInWorld(origin_in_world),
    _timer(std::make_unique<G4Timer>()),

    _processInfo(phys_process_info),
    _g4InternalFiltering(conf.G4InteralFiltering())
  {
    auto c = global_c;
    if(conf.SimParticlePrinter(c)) {
      simParticlePrinter_ = SimParticleCollectionPrinter(c);
    }

    //G4SDManager* SDman = G4SDManager::GetSDMpointer();
    //SDman->ListTree();
  }

  Mu2eG4EventAction::~Mu2eG4EventAction()
  {}


  void Mu2eG4EventAction::BeginOfEventAction(const G4Event *evt)
  {
    art::Event *_artEvent = perThreadObjects_->artEvent;
    SimParticleHelper *_spHelper = &perThreadObjects_->simParticleHelper;
    SimParticlePrimaryHelper *_parentHelper = &perThreadObjects_->simParticlePrimaryHelper;

    // local Mu2e timer, almost equal to time of G4EventManager::ProcessOneEvent()
    _timer->Start();

    art::Handle<SimParticleCollection> inputSimHandle;
    if(art::InputTag() != multiStagePars_.inputSimParticles()) {
      _artEvent->getByLabel(multiStagePars_.inputSimParticles(), inputSimHandle);

      if(!inputSimHandle.isValid()) {
        throw cet::exception("CONFIG")
          << "Error retrieving inputSimParticles for "
          << multiStagePars_.inputSimParticles() <<"\n";
      }
    }

    art::Handle<MCTrajectoryCollection> inputMCTrajectoryHandle;
    if(art::InputTag() != multiStagePars_.inputMCTrajectories()) {
      _artEvent->getByLabel(multiStagePars_.inputMCTrajectories(), inputMCTrajectoryHandle);

      if(!inputMCTrajectoryHandle.isValid()) {
        throw cet::exception("CONFIG")
          << "Error retrieving inputMCTrajectories for "
          << multiStagePars_.inputMCTrajectories() <<"\n";
      }
    }

    //these are OK, nothing put into or defined for event
    _sensitiveDetectorHelper->createProducts(*_artEvent, *_spHelper);

    //note: these calls to finishConstruction() must happen AFTER the call to userDetector->Construct().
    //This call occurs in G4RunManager::InitializeGeometry(), which occurs in _runManager->Initialize() in Mu2eG4_module.
    //We cannot put these calls to finishConstruction() in ActionInitialization::Build(), because in sequential mode
    //Build() is called at the call to _runManager->SetUserInitialization(actioninit), which happens BEFORE _runManager->Initialize().

    _stackingCuts->finishConstruction(_originInWorld);
    _steppingCuts->finishConstruction(_originInWorld);
    _commonCuts->finishConstruction(_originInWorld);

    _stackingCuts->beginEvent(*_artEvent, *_spHelper);
    _steppingCuts->beginEvent(*_artEvent, *_spHelper);
    _commonCuts->beginEvent(*_artEvent, *_spHelper);

    _trackingAction->beginEvent(inputSimHandle, inputMCTrajectoryHandle, *_spHelper,
                                *_parentHelper, *perThreadObjects_->mcTrajectories,
                                *perThreadObjects_->simRemapping);

    _steppingAction->BeginOfEvent(*perThreadObjects_->tvd_collection, *_spHelper);

    _sensitiveDetectorHelper->updateSensitiveDetectors(*_processInfo, *_spHelper);

    if (_sensitiveDetectorHelper->getExtMonFNALPixelSD()) {
      perThreadObjects_->extMonFNALHits = unique_ptr<ExtMonFNALSimHitCollection>( new ExtMonFNALSimHitCollection );
      _sensitiveDetectorHelper->getExtMonFNALPixelSD()->beforeG4Event(perThreadObjects_->extMonFNALHits.get(), *_spHelper);
    }

  }


  void Mu2eG4EventAction::EndOfEventAction(const G4Event *evt)
  {
    // Run self consistency checks if enabled.
    _trackingAction->endEvent(*perThreadObjects_->simPartCollection);

    _timer->Stop();

    simParticlePrinter_.print(std::cout, *perThreadObjects_->simPartCollection);

    // Pass data products to the module to put into the event
    bool event_passes = false;

    //if internal G4 filtering is turned OFF, event passes
    //otherwise event must pass the StepPoints momentum and Tracker StepPoints filters
    if (_g4InternalFiltering) {
      if (_sensitiveDetectorHelper->filterStepPointMomentum() && _sensitiveDetectorHelper->filterTrackerStepPoints()) {
        event_passes = true;
      }
    } else {
      event_passes = true;
    }


    if (event_passes) {

      // Populate the output data products.
      // Fill the status object.
      float cpuTime  = _timer->GetSystemElapsed() + _timer->GetUserElapsed();

      int status(0);
      if ( _steppingAction->nKilledStepLimit() > 0 ||
           _trackingAction->nKilledByFieldPropagator() > 0 ) {
        status =  1;
      }
      if ( _trackingAction->overflowSimParticles() ) status = 10;

      perThreadObjects_->statG4 = std::make_unique<StatusG4>(status,
                                                             _trackingAction->nG4Tracks(),
                                                             _trackingAction->overflowSimParticles(),
                                                             _steppingAction->nKilledStepLimit(),
                                                             _trackingAction->nKilledByFieldPropagator(),
                                                             cpuTime,
                                                             _timer->GetRealElapsed()
                                                             );

      _sensitiveDetectorHelper->insertSDDataIntoPerThreadStorage(perThreadObjects_);
      _stackingCuts->insertCutsDataIntoPerThreadStorage(perThreadObjects_);
      _steppingCuts->insertCutsDataIntoPerThreadStorage(perThreadObjects_);
      _commonCuts->insertCutsDataIntoPerThreadStorage(perThreadObjects_);

    }//event_passes cuts
    else {
      //there is no need to clear SD data here, since it is done in the call to
      //_sensitiveDetectorHelper->createProducts in the BeginOfEventAction above

      _stackingCuts->deleteCutsData();
      _steppingCuts->deleteCutsData();
      _commonCuts->deleteCutsData();
    }//else put NULL ptrs into the stash

  }//EndOfEventAction

} // end namespace mu2e
