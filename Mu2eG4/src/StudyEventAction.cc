//
// G4 begin and end of event actions for Mu2e test environment
//
// $Id: StudyEventAction.cc,v 1.1 2012/11/16 23:53:04 genser Exp $
// $Author: genser $
// $Date: 2012/11/16 23:53:04 $
//
// Original author KLG
//

#include "Mu2eG4/inc/StudyEventAction.hh"
#include "Mu2eG4/inc/StudySteppingAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

namespace mu2e {

StudyEventAction::StudyEventAction(StudySteppingAction *stepping_action) {
  _stepping = stepping_action;
}

StudyEventAction::~StudyEventAction()
{}


void StudyEventAction::BeginOfEventAction(const G4Event*) {
  //  _stepping->BeginOfEvent();
}


void StudyEventAction::EndOfEventAction(const G4Event* evt)
{
  _stepping->EndOfEvent();

  // Change  G4_plugin so that we copy to the art::event from here.

}

} // end namespace mu2e
