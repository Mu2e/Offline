#ifndef Mu2eG4_EventAction_hh
#define Mu2eG4_EventAction_hh
//
// G4 begin and end of event actions for Mu2e.
//
// $Id: EventAction.hh,v 1.4 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
// Original author Rob Kutschke
//

#include "G4UserEventAction.hh"


class G4Event;

namespace mu2e {

class SteppingAction;

class EventAction : public G4UserEventAction
{
  public:
    EventAction(SteppingAction*);
   ~EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

  private:
    SteppingAction * _stepping;
};

}  // end namespace mu2e
#endif /* Mu2eG4_EventAction_hh */


