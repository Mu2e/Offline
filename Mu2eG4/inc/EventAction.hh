#ifndef EventAction_h
#define EventAction_h 1
//
// G4 begin and end of event actions for Mu2e.
//
// $Id: EventAction.hh,v 1.2 2010/09/20 02:57:05 logash Exp $
// $Author: logash $ 
// $Date: 2010/09/20 02:57:05 $
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
#endif

    
