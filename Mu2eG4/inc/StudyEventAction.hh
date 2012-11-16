#ifndef Mu2eG4_StudyEventAction_hh
#define Mu2eG4_StudyEventAction_hh
//
// G4 begin and end of event actions for Mu2e test environment
//
// $Id: StudyEventAction.hh,v 1.1 2012/11/16 23:29:03 genser Exp $
// $Author: genser $
// $Date: 2012/11/16 23:29:03 $
//
// Original author KLG
//

#include "G4UserEventAction.hh"

class G4Event;

namespace mu2e {

class StudySteppingAction;

class StudyEventAction : public G4UserEventAction
{
  public:
    StudyEventAction(StudySteppingAction*);
   ~StudyEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

  private:
    StudySteppingAction * _stepping;
};

}  // end namespace mu2e
#endif /* Mu2eG4_StudyEventAction_hh */


