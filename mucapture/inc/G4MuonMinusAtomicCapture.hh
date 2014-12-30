#ifndef G4MuMinusAtomicCapture_h
#define G4MuMinusAtomicCapture_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4VRestProcess.hh"
#include "G4ProcessType.hh"
#include "G4ParticleDefinition.hh"
#include "G4VParticleChange.hh"
#include "G4ElementSelector.hh"

class G4MuonMinusAtomicCapture : public G4VRestProcess {

public:
  G4MuonMinusAtomicCapture(G4String const& name = "muMinusAtomicCapture",
			   G4ProcessType type = fElectromagnetic);							    // type? 

  virtual ~G4MuonMinusAtomicCapture();
  
  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  virtual G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&); 
  
  G4double GetMeanLifeTime(const G4Track&, G4ForceCondition*){return 0;};



  void PreparePhysicsTable(const G4ParticleDefinition& p);

  void BuildPhysicsTable(const G4ParticleDefinition& p);




private:
  G4ElementSelector *pSelector;
};


#endif
