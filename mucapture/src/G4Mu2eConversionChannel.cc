// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, March 22 include
// ----------------------------------------------------------------

#include "G4Mu2eConversionChannel.hh"

#include "G4Electron.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include <sstream>

G4Mu2eConversionChannel::
G4Mu2eConversionChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "Mu2e Conversion", captureRate, verboseLevel),
  psdc(GetPSDC(p)) {
}

G4Mu2eConversionChannel::~G4Mu2eConversionChannel(){
  delete psdc;
}

G4DecayProducts* G4Mu2eConversionChannel::CaptureIt(G4DynamicParticle const* pParticle){
  G4DecayProducts *products = psdc->DecayIt(pParticle->GetMass());
  products->SetParentParticle(*pParticle);
  return products;
}

G4PhaseSpaceDecayChannel* G4Mu2eConversionChannel::
GetPSDC(G4MuAtom const* pParticle){
  G4int const Z = pParticle->GetAtomicNumber();
  G4int const A = pParticle->GetAtomicMass();
  G4IonTable *table = G4ParticleTable::GetParticleTable()->GetIonTable();
  // make sure the ion is created
  G4ParticleDefinition *ion = table->GetIon(Z,A);
  G4String const IonName = ion->GetParticleName();
  return new 
    G4PhaseSpaceDecayChannel(pParticle->GetParticleName(),
			     1.,2,"e-",IonName.c_str());
}
