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
    trajectoryControl_(conf.TrajectoryControl()),
    simParticlePrinter_(),
    timeVDtimes_(conf.SDConfig().TimeVD().times()),
    multiStagePars_(conf),
    _trackingAction(tracking_action),
    _steppingAction(stepping_action),
    _sensitiveDetectorHelper(sensitive_detectorhelper),
    _stackingCuts(&stacking_cuts),
    _steppingCuts(&stepping_cuts),
    _commonCuts(&common_cuts),
    _originInWorld(origin_in_world),
    _tvdOutputName(StepInstanceName::timeVD),
    _timer(std::make_unique<G4Timer>()),

    _spHelper(),
    _parentHelper(),
    _processInfo(phys_process_info),
    _artEvent(),
    _g4InternalFiltering(conf.G4InteralFiltering())
  {
    SimParticleCollectionPrinter::Config c;
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
    setEventData();

    _spHelper = perThreadObjects_->simParticleHelper;
    _parentHelper = perThreadObjects_->simParticlePrimaryHelper;

    // local Mu2e timer, almost equal to time of G4EventManager::ProcessOneEvent()
    _timer->Start();

    simParticles = unique_ptr<SimParticleCollection>( new SimParticleCollection );
    tvdHits = unique_ptr<StepPointMCCollection>( new StepPointMCCollection );
    mcTrajectories = unique_ptr<MCTrajectoryCollection>( new MCTrajectoryCollection );
    simsRemap = unique_ptr<SimParticleRemapping>( new SimParticleRemapping );
    extMonFNALHits = unique_ptr<ExtMonFNALSimHitCollection>( new ExtMonFNALSimHitCollection );

    //these will NEVER be run in MT mode, so we don't need a mutex and lock
    //on the _artEvent->getByLabel call, even though the call is not thread-safe
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
                                *_parentHelper, *mcTrajectories, *simsRemap);

    _steppingAction->BeginOfEvent(*tvdHits, *_spHelper);

    _sensitiveDetectorHelper->updateSensitiveDetectors(*_processInfo, *_spHelper);

    if (_sensitiveDetectorHelper->getExtMonFNALPixelSD()) {
      _sensitiveDetectorHelper->getExtMonFNALPixelSD()->beforeG4Event(extMonFNALHits.get(), *_spHelper);
    }

  }


  void Mu2eG4EventAction::EndOfEventAction(const G4Event *evt)
  {
    // Run self consistency checks if enabled.
    _trackingAction->endEvent(*simParticles);

    _timer->Stop();

    // Populate the output data products.
    // Fill the status object.
    float cpuTime  = _timer->GetSystemElapsed() + _timer->GetUserElapsed();

    int status(0);
    if ( _steppingAction->nKilledStepLimit() > 0 ||
         _trackingAction->nKilledByFieldPropagator() > 0 ) {
      status =  1;
    }
    if ( _trackingAction->overflowSimParticles() ) status = 10;

    auto g4stat = std::make_unique<StatusG4>(status,
                                             _trackingAction->nG4Tracks(),
                                             _trackingAction->overflowSimParticles(),
                                             _steppingAction->nKilledStepLimit(),
                                             _trackingAction->nKilledByFieldPropagator(),
                                             cpuTime,
                                             _timer->GetRealElapsed()
                                             );

    simParticlePrinter_.print(std::cout, *simParticles);

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
      perThreadObjects_->insertSimsAndStatusData(std::move(g4stat), std::move(simParticles));

      _sensitiveDetectorHelper->insertSDDataIntoPerThreadStorage(perThreadObjects_);
      _stackingCuts->insertCutsDataIntoPerThreadStorage(perThreadObjects_);
      _steppingCuts->insertCutsDataIntoPerThreadStorage(perThreadObjects_);
      _commonCuts->insertCutsDataIntoPerThreadStorage(perThreadObjects_);

      if(!timeVDtimes_.empty()) {
        perThreadObjects_->insertTVDHits(std::move(tvdHits), _tvdOutputName.name());
      }

      if(trajectoryControl_.produce()) {
        perThreadObjects_->insertMCTrajectoryCollection(std::move(mcTrajectories));
      }

      if(multiStagePars_.multiStage()) {
        perThreadObjects_->insertSimsRemapping(std::move(simsRemap));
      }

      if(_sensitiveDetectorHelper->extMonPixelsEnabled()) {
        perThreadObjects_->insertExtMonFNALSimHits(std::move(extMonFNALHits));
      }

    }//event_passes cuts
    else {
      simParticles = nullptr;
      tvdHits = nullptr;
      mcTrajectories = nullptr;
      simsRemap = nullptr;
      extMonFNALHits = nullptr;

      //there is no need to clear SD data here, since it is done in the call to
      //_sensitiveDetectorHelper->createProducts in the BeginOfEventAction above

      _stackingCuts->deleteCutsData();
      _steppingCuts->deleteCutsData();
      _commonCuts->deleteCutsData();
    }//else put NULL ptrs into the stash

  }//EndOfEventAction


  void Mu2eG4EventAction::setEventData()
  {
    _artEvent = perThreadObjects_->artEvent;
  }

} // end namespace mu2e
