#ifndef WLSRunAction_h
#define WLSRunAction_h 1

#include "globals.hh"

#include "G4UserRunAction.hh"

class G4Run;

class WLSRunAction : public G4UserRunAction
{
  public:

    WLSRunAction(int mode);
    ~WLSRunAction();

  public:

    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

    void  SetRndmFreq(G4int val) { saveRndm = val; }
    G4int GetRndmFreq()          { return saveRndm; }

    inline void SetAutoSeed (const G4bool val) { autoSeed = val; }

  private:
 
    G4int saveRndm;
    G4bool autoSeed;

    int _mode;

};

#endif
