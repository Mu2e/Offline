//
// G4 begin and end of event actions for Mu2e.
//
// $Id: EventAction.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

namespace mu2e {

EventAction::EventAction()
{}

EventAction::~EventAction()
{}

 
void EventAction::BeginOfEventAction(const G4Event*)
{}

 
void EventAction::EndOfEventAction(const G4Event* evt)
{

  // Change  G4_plugin so that we copy to the edm::event from here.

}

} // end namespace mu2e
