#ifndef G4GenericMuAtom_h
#define G4GenericMuAtom_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "G4MuAtom.hh"

class G4GenericMuAtom : public G4MuAtom
{
public:
  static G4GenericMuAtom* Definition();
  static G4GenericMuAtom* GenericMuAtom();
  static G4GenericMuAtom* GenericMuAtomDefinition();

private:
  static G4GenericMuAtom* theInstance;
  G4GenericMuAtom(const G4String&  aName,  
		  G4double         mass,     
		  G4int            Z,
		  G4int            A, 
		  G4int            iSpin,
		  G4int            encoding);

  ~G4GenericMuAtom(){}
};

#endif
