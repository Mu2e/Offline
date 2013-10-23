// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

#include "G4GenericMuAtom.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"

G4GenericMuAtom* G4GenericMuAtom::theInstance = 0;


G4GenericMuAtom* G4GenericMuAtom::Definition(){
  return GenericMuAtomDefinition();
}

G4GenericMuAtom* G4GenericMuAtom::GenericMuAtom(){
  return GenericMuAtomDefinition();
}

G4GenericMuAtom* G4GenericMuAtom::GenericMuAtomDefinition(){
  if( theInstance != 0 )
    return theInstance;

  G4String const name("GenericMuAtom");
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  // FIXME ... should this be a dynamic_cast or a static_cast
  G4GenericMuAtom* anInstance = static_cast<G4GenericMuAtom*>(pTable->FindParticle(name));
  if( anInstance == 0 )
    theInstance = new G4GenericMuAtom(name, 10.*CLHEP::GeV, 1, 1, 0, 0);
  return theInstance;
}


G4GenericMuAtom::G4GenericMuAtom(const G4String&  aName,  
				 G4double         mass,     
				 G4int            Z,
				 G4int            A, 
				 G4int            iSpin,
				 G4int            encoding) :
  G4MuAtom(aName, mass, Z, A, iSpin, encoding){}
