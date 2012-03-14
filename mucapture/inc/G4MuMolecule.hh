#ifndef G4MuMolecule_h
#define G4MuMolecule_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, March 23 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4ParticleDefinition.hh"

#include "G4MuMoleculeCaptureKineticsTable.hh"

class G4MuMolecule : public G4ParticleDefinition {
  // Class Description:

  // This class declares a new particle type for describing light
  // molecules bound by negative muons ("muonic molecules").  In
  // addition to the G4DecayTable, this class adds a
  // G4MuMoleculeKineticsTable with similar semantics to the G4MuAtom
  // G4MuAtomCaptureKineticsTable 

public:

  virtual ~G4MuMolecule();

  // G4MuMolecules must be explicitly intantiated
  static G4MuMolecule* MuMolecule(G4int Z1, G4int A1, G4int Z2, G4int A2, G4int iSpin=0);
  static G4MuMolecule* 
  MuMoleculeDefinition(G4int Z1, G4int A1, G4int Z2, G4int A2, G4int iSpin=0);
  static G4MuMolecule* Definition(G4int Z1, G4int A1, G4int Z2, G4int A2, G4int iSpin=0);

  // get/set decay table ... we can use the normal one here

  G4MuMoleculeCaptureKineticsTable*
  CaptureKineticsTable(G4MuMoleculeCaptureKineticsTable*);
  G4MuMoleculeCaptureKineticsTable*
  CaptureKineticsTable() const { return fCaptureTable; }

  // get the capture rate ... from kinetics table
  G4double GetCaptureRate() const;

  // IMPLEMENT_ME The design calls for state and species selection
  // models ... we may also want a charge model and a lifetime model.
  // Perhaps we need a G4MuMoleculeTable, or perhaps the models just
  // get attached to the MuMolecule directly?  there aren't more than
  // a handful of molecules to worry about...
  // Furthermore, is there a cascade to worry about?

protected:
  G4MuMolecule(G4String const& name,
	       G4double mass,
	       G4double Z1, G4double A1,
	       G4double Z2, G4double A2,
	       G4int iSpin, G4int encoding);

private:
  
  // order these such that A1<=A2
  G4int Z1, A1, Z2, A2;

  G4MuMoleculeCaptureKineticsTable* fCaptureTable;

  static void OrderParams(G4int& Z1, G4int& A1, G4int& Z2, G4int& A2);
  static G4String MakeName(G4int Z1, G4int A1, G4int Z2, G4int A2, G4int iSpin);
  static G4String MakeName(G4int Z, G4int A);
  static G4int MakeEncoding(G4int Z1, G4int A1, G4int Z2, G4int A2);

};

#endif
