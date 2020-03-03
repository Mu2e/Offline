//
// ActionInitialization.cc provides implementation of Mu2e G4's built-in action initialization.
//
// Author: Lisa Goodenough
// Date: 2017/05/08
//
//

//Mu2e includes
#include "Mu2eG4/inc/ActionInitialization.hh"
#include "Mu2eG4/inc/PrimaryGeneratorAction.hh"
#include "Mu2eG4/inc/Mu2eG4StackingAction.hh"
#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/Mu2eG4SteppingAction.hh"
#include "Mu2eG4/inc/Mu2eG4EventAction.hh"
#include "Mu2eG4/inc/Mu2eG4RunAction.hh"
#include "Mu2eG4/inc/Mu2eG4MasterRunAction.hh"
#include "Mu2eG4/inc/ExtMonFNALPixelSD.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "Mu2eG4/inc/IMu2eG4Cut.hh"
#include "Mu2eG4/inc/Mu2eG4PerThreadStorage.hh"
#include "Mu2eG4/inc/SteppingVerbose.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"

//C++ includes
#include <iostream>

namespace mu2e {

  ActionInitialization::ActionInitialization(const Mu2eG4Config::Top& conf,
                                             SensitiveDetectorHelper* sensitive_detectorhelper,
                                             Mu2eG4PerThreadStorage* per_thread_storage,
                                             PhysicalVolumeHelper* phys_volume_helper,
                                             CLHEP::Hep3Vector const& origin_in_world,
                                             unsigned stage_offset_for_tracking_action
                                             )
  :
  G4VUserActionInitialization(),
    conf_(conf),
    trajectoryControl_(conf.TrajectoryControl()),
    timeVDtimes_(conf.SDConfig().TimeVD().times()),
    mu2eLimits_(conf.ResourceLimits()),

    stackingCuts_(createMu2eG4Cuts(conf.Mu2eG4StackingOnlyCut.get<fhicl::ParameterSet>(), mu2eLimits_)),
    steppingCuts_(createMu2eG4Cuts(conf.Mu2eG4SteppingOnlyCut.get<fhicl::ParameterSet>(), mu2eLimits_)),
    commonCuts_(createMu2eG4Cuts(conf.Mu2eG4CommonCut.get<fhicl::ParameterSet>(), mu2eLimits_)),

    sensitiveDetectorHelper_(sensitive_detectorhelper),
    perThreadStorage_(per_thread_storage),
    physVolHelper_(phys_volume_helper),
    physicsProcessInfo_(),

    originInWorld_(origin_in_world),
    stageOffset_(stage_offset_for_tracking_action)
  {}

  ActionInitialization::~ActionInitialization()
  {}

  //nothing to do, this is only for the Master Thread
  void ActionInitialization::BuildForMaster() const
  {}


  // used for defining user action classes in sequential mode.
  void ActionInitialization::Build() const
  {
    IMu2eG4Cut& stacking_Cuts = *stackingCuts_.get();
    IMu2eG4Cut& stepping_Cuts = *steppingCuts_.get();
    IMu2eG4Cut& common_Cuts = *commonCuts_.get();

    PrimaryGeneratorAction* genAction = new PrimaryGeneratorAction(conf_.debug(), perThreadStorage_);
    SetUserAction(genAction);

    Mu2eG4SteppingAction* steppingAction = new Mu2eG4SteppingAction(conf_.debug(),
                                                                    timeVDtimes_,
                                                                    stepping_Cuts,
                                                                    common_Cuts,
                                                                    trajectoryControl_,
                                                                    mu2eLimits_);
    SetUserAction(steppingAction);

    SetUserAction( new Mu2eG4StackingAction(stacking_Cuts, common_Cuts) );

    TrackingAction* trackingAction = new TrackingAction(conf_,
                                                        steppingAction,
                                                        stageOffset_,
                                                        trajectoryControl_,
                                                        mu2eLimits_);
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
                                         stacking_Cuts,
                                         stepping_Cuts,
                                         common_Cuts,
                                         perThreadStorage_,
                                         &physicsProcessInfo_,
                                         originInWorld_) );

  }//Build()


  G4VSteppingVerbose* ActionInitialization::InitializeSteppingVerbose() const
  {
    return new SteppingVerbose;
  }


} // end namespace mu2e
