#ifndef WLSRunAction_h
#define WLSRunAction_h 1

#include "globals.hh"

#include "G4UserRunAction.hh"

class G4Run;

class WLSRunAction : public G4UserRunAction
{
  public:

    WLSRunAction();
    ~WLSRunAction();

    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

};

#endif
