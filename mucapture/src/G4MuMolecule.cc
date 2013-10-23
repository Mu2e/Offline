// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, March 23 2010
// ----------------------------------------------------------------

#include "G4MuMolecule.hh"
#include "G4MuAtom.hh"
#include "G4Ions.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4IonTable.hh"

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

#include <sstream>
#include <utility>

G4MuMolecule::~G4MuMolecule(){

}

G4MuMolecule* 
G4MuMolecule::MuMolecule(G4int Z1, G4int A1, G4int Z2, G4int A2, G4int iSpin){
  return Definition(Z1,A1,Z2,A2,iSpin);
}

G4MuMolecule* 
G4MuMolecule::MuMoleculeDefinition(G4int Z1, G4int A1, G4int Z2, G4int A2, G4int iSpin){ 
  return Definition(Z1,A1,Z2,A2,iSpin);
}

G4MuMolecule* 
G4MuMolecule::Definition(G4int Z1, G4int A1, G4int Z2, G4int A2, G4int iSpin){
  OrderParams(Z1,A1,Z2,A2);
  G4ParticleTable *table = G4ParticleTable::GetParticleTable();
  G4String const name = G4MuMolecule::MakeName(Z1,A1,Z2,A2,iSpin);
  G4MuMolecule *mumol = static_cast<G4MuMolecule*>( table->FindParticle(name) );
  if( mumol == 0 ){ // need to create it
    G4String const name = G4MuMolecule::MakeName(Z1,A1,Z2,A2,iSpin);
    G4int const encoding = G4MuMolecule::MakeEncoding(Z1,A1,Z2,A2);
    // FIXME ... ought to be able to get this from somewhere!
    G4double const muon_mass_c2 = 0.1056584*CLHEP::GeV;
    G4double mass = G4NucleiProperties::GetNuclearMass(A1,Z1) + 
      G4NucleiProperties::GetNuclearMass(A2,Z2) +
      muon_mass_c2;
    // FIXME ... binding energy?
    mumol = new G4MuMolecule(name, mass, Z1,A1,Z2,A2,iSpin,encoding); // leak
  }
  return mumol;
}

// FIXME ... lots of things about this are wrong ... the charge, for
// starters ...
G4MuMolecule::G4MuMolecule(G4String const& name,
			   G4double mass,
			   G4double Z1, G4double A1,
			   G4double Z2, G4double A2,
			   G4int iSpin, G4int encoding) :
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding
  G4ParticleDefinition( name,mass,0,(Z1+Z2-1)*CLHEP::eplus,iSpin,0,0,0,0,0,
			"MuMolecule",-1,A1+A2,encoding,false,0.,0,
			false,"",0),
  Z1(Z1), A1(A1), Z2(Z2), A2(A2), fCaptureTable(0) {}

G4MuMoleculeCaptureKineticsTable*
G4MuMolecule::CaptureKineticsTable(G4MuMoleculeCaptureKineticsTable* table){
  std::swap(table, fCaptureTable);
  return table;
}

G4double G4MuMolecule::GetCaptureRate() const {
  if( !fCaptureTable )
    return 0;
  return fCaptureTable->GetTotalRate();
}

// for reasons that are clear only to me in a former time, we want the
// lighter one on the left of the name
void G4MuMolecule::OrderParams(G4int& Z1, G4int& A1, G4int& Z2, G4int& A2){
  if(A1>A2){
    std::swap(A1,A2);
    std::swap(Z1,Z2);
  }
}

G4String G4MuMolecule::MakeName(G4int Z1, G4int A1, G4int Z2, G4int A2, G4int iSpin){ 
  // nucleus_mu_nucleus_(iSpin/2)
  G4String i1 = G4MuMolecule::MakeName(Z1,A1), i2 = G4MuMolecule::MakeName(Z2,A2);
  std::ostringstream o;
  o << i1 << "_mu_" << i2;
  if( iSpin!=0 ){
    if( iSpin%2 == 0 )
      o << '_' << iSpin/2;
    else
      o << '_' << iSpin << "/2";
  }
  return o.str();
}

G4String G4MuMolecule::MakeName(G4int Z, G4int A){
  G4String name;
  if( Z==1 ){
    if( A == 1 )
      name = "P";
    else if ( A == 2 )
      name = "D";
    else if ( A == 3 )
      name = "T";
  } else {
    std::ostringstream o;
    o << G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonName(Z,A,0.0);
    name = o.str();
    G4int pos = name.find('[');
    // string [nnnn]
    name = name.substr(0,pos);
  }
  return name;
}

// FIXME  ... I guess I have to invent a good encoding
G4int G4MuMolecule::MakeEncoding(G4int Z1, G4int A1, G4int /*Z2*/, G4int /*A2*/){ 
  G4int encoding = 
    G4ParticleTable::GetParticleTable()->GetIonTable()->GetNucleusEncoding(Z1, A1, 0);
  encoding += 2000000000;
  return encoding;
}
