// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, March 22 include
// ----------------------------------------------------------------

#include "G4Mu2eOrbitalConversionChannel.hh"

#include "G4MuAtom.hh"
#include "G4Electron.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

// base class

G4Mu2eOrbitalConversionChannel::
G4Mu2eOrbitalConversionChannel(G4MuAtom const* p,
			       G4double captureRate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "Mu2e Orbital Conversion", 
				  captureRate, verboseLevel),
  psdc( GetPSDC(p) )
{}


// FIXME/IMPLEMENT_ME ... this is certainly wrong :-) To do this
// right, I need to account for the orbital velocities of the mu- and
// e-, though their binding energies, get their CM frame, do a boost
// .... look at G4MuAtomDIOChannel, I believe (after confirming that
// that works properly :-)
G4DecayProducts* 
G4Mu2eOrbitalConversionChannel::CaptureIt(G4DynamicParticle const* pParticle){
  G4DecayProducts *products =
    psdc->DecayIt(pParticle->GetDefinition()->GetPDGMass());
  products->SetParentParticle(*pParticle);
  return products;
}

G4PhaseSpaceDecayChannel* G4Mu2eOrbitalConversionChannel::
GetPSDC(G4MuAtom const* pParticle){
  G4int const Z = pParticle->GetAtomicNumber();
  G4int const A = pParticle->GetAtomicMass();
  G4IonTable *table = G4ParticleTable::GetParticleTable()->GetIonTable();
  // force Ion creation
  G4ParticleDefinition *ion = table->GetIon(Z,A);
  G4String const IonName = ion->GetParticleName();

  return
    new G4PhaseSpaceDecayChannel(pParticle->GetParticleName(),
				 1.,3,"e-", "e-", IonName.c_str());
}
