#ifndef EventAction_h
#define EventAction_h 1
//
// G4 begin and end of event actions for Mu2e.
//
// $Id: EventAction.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include "G4UserEventAction.hh"


class G4Event;

namespace mu2e {

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
   ~EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
};

}  // end namespace mu2e
#endif

    
