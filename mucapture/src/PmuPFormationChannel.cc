// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, March 12 2010
// ----------------------------------------------------------------

#include "PmuPFormationChannel.hh"

#include "G4MuMolecule.hh"
#include "G4ThreeVector.hh"

#include "RandomUtilities.hh"

#if G4VERSION<4095
#include <strstream>
#endif

// base class

PmuPFormationChannel::
PmuPFormationChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "P_mu_P formation", captureRate, verboseLevel) {
  CheckIsApplicable();
}

G4DecayProducts* PmuPFormationChannel::CaptureIt(G4DynamicParticle const* pParticle){
  G4DecayProducts *products = new G4DecayProducts;
  // Set the parent
  products->SetParentParticle(*pParticle);

  G4DynamicParticle *pmup = 
    new G4DynamicParticle(G4MuMolecule::Definition(1,1,1,1), G4ThreeVector());
  products->PushProducts(pmup);

  return products;
}


void PmuPFormationChannel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,4) != "mu_P" ){
#if G4VERSION<4095
    std::ostrstream ed;
#else
    G4ExceptionDescription ed;
#endif
    ed << "Channel " << GetChannelName() << " is not applicable to particle "
       << part->GetParticleName();
    G4Exception("PmuPFormationChannel::CheckIsApplicable",
                "PMUPF0001", FatalException, 
#if G4VERSION<4095
                ed.str());
#else
                ed);
#endif
  }
}
