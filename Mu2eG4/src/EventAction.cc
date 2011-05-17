//
// G4 begin and end of event actions for Mu2e.
//
// $Id: EventAction.cc,v 1.3 2011/05/17 15:36:00 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:36:00 $
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/EventAction.hh"
#include "Mu2eG4/inc/SteppingAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

namespace mu2e {

EventAction::EventAction(SteppingAction *stepping_action) { 
  _stepping = stepping_action;
}

EventAction::~EventAction()
{}

 
void EventAction::BeginOfEventAction(const G4Event*) {
  _stepping->BeginOfEvent();
}

 
void EventAction::EndOfEventAction(const G4Event* evt)
{
  _stepping->EndOfEvent();

  // Change  G4_plugin so that we copy to the art::event from here.

}

} // end namespace mu2e
