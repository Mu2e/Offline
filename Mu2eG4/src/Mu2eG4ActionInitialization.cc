//
// Mu2eG4ActionInitialization.cc provides implementation of Mu2e G4's built-in action initialization.
//
// Author: Lisa Goodenough
// Date: 2017/05/08
//
//

//Mu2e includes
#include "Mu2eG4/inc/Mu2eG4ActionInitialization.hh"
#include "Mu2eG4/inc/Mu2eG4PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/Mu2eG4StackingAction.hh"
#include "Mu2eG4/inc/Mu2eG4TrackingAction.hh"
#include "Mu2eG4/inc/Mu2eG4SteppingAction.hh"
#include "Mu2eG4/inc/Mu2eG4EventAction.hh"
#include "Mu2eG4/inc/Mu2eG4RunAction.hh"
#include "Mu2eG4/inc/Mu2eG4MasterRunAction.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/IMu2eG4Cut.hh"
#include "Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Mu2eG4/inc/Mu2eG4SteppingVerbose.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"

//C++ includes
#include <iostream>

namespace mu2e {

  Mu2eG4ActionInitialization::Mu2eG4ActionInitialization(const Mu2eG4Config::Top& conf,
                                             SensitiveDetectorHelper* sensitive_detectorhelper,
                                             Mu2eG4PerThreadStorage* per_thread_storage,
                                             PhysicalVolumeHelper* phys_volume_helper,
                                             CLHEP::Hep3Vector const& origin_in_world
                                             )
  :
  G4VUserActionInitialization(),
    conf_(conf),
    timeVDtimes_(conf.SDConfig().TimeVD().times()),

    sensitiveDetectorHelper_(sensitive_detectorhelper),
    perThreadStorage_(per_thread_storage),
    physVolHelper_(phys_volume_helper),
    physicsProcessInfo_(),

    originInWorld_(origin_in_world)
  {}

  Mu2eG4ActionInitialization::~Mu2eG4ActionInitialization()
  {}

  //nothing to do, this is only for the Master Thread
  void Mu2eG4ActionInitialization::BuildForMaster() const
  {}


  // used for defining user action classes in sequential mode.
  void Mu2eG4ActionInitialization::Build() const
  {
    Mu2eG4PrimaryGeneratorAction* genAction = new Mu2eG4PrimaryGeneratorAction(conf_.debug(), perThreadStorage_);
    SetUserAction(genAction);

    Mu2eG4SteppingAction* steppingAction = new Mu2eG4SteppingAction(conf_.debug(),
                                                                    timeVDtimes_,
                                                                    *perThreadStorage_->steppingCuts,
                                                                    *perThreadStorage_->commonCuts,
                                                                    perThreadStorage_->ioconf.trajectoryControl(),
                                                                    perThreadStorage_->ioconf.mu2elimits());
    SetUserAction(steppingAction);

    SetUserAction( new Mu2eG4StackingAction(*perThreadStorage_->stackingCuts, *perThreadStorage_->commonCuts) );

    Mu2eG4TrackingAction* trackingAction = new Mu2eG4TrackingAction(conf_,
                                                        steppingAction,
                                                        perThreadStorage_);
    SetUserAction(trackingAction);

    SetUserAction( new Mu2eG4RunAction(conf_.debug(),
                                       originInWorld_,
                                       physVolHelper_,
                                       &physicsProcessInfo_,
                                       trackingAction,
                                       steppingAction,
                                       sensitiveDetectorHelper_) );


    SetUserAction( new Mu2eG4EventAction(conf_,
                                         trackingAction,
                                         steppingAction,
                                         sensitiveDetectorHelper_,
                                         perThreadStorage_,
                                         &physicsProcessInfo_,
                                         originInWorld_) );

  }//Build()


  G4VSteppingVerbose* Mu2eG4ActionInitialization::InitializeSteppingVerbose() const
  {
    return new Mu2eG4SteppingVerbose;
  }


} // end namespace mu2e
