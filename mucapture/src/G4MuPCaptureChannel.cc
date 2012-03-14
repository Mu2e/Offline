// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, March 12 2010
// ----------------------------------------------------------------

#include "G4MuPCaptureChannel.hh"

#include "G4Neutron.hh"
#include "G4NeutrinoMu.hh"
#include "G4ThreeVector.hh"
#include "G4MuAtomTable.hh"

#include "RandomUtilities.hh"

#include <strstream>

// base class

G4MuPCaptureChannel::
G4MuPCaptureChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_P capture", captureRate, verboseLevel) {
  CheckIsApplicable();
}

G4MuPCaptureChannel::
G4MuPCaptureChannel(G4MuAtom const* p, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_P capture", 0., verboseLevel) {
  CheckIsApplicable();
  G4int const Z = p->GetAtomicNumber();
  G4int const A = p->GetAtomicMass();
  G4int const iSpin = p->GetPDGiSpin();
  G4MuAtomTable const* table = G4MuAtomTable::GetInstance();
  G4MuAtomCaptureRateModel const* model = table->CaptureRateModel(Z,A);
  rate = model->GetCaptureRate(Z,A,iSpin);
}

// FIXME ... Reimplement with psdc?
G4DecayProducts* G4MuPCaptureChannel::CaptureIt(G4DynamicParticle const* pParticle){
  G4DecayProducts *products = new G4DecayProducts;

  // Set the parent
  products->SetParentParticle(*pParticle);

  // start calculating  
  G4double const pmass = pParticle->GetMass();
  G4double const nmass = G4Neutron::Definition()->GetPDGMass();

  G4double const nuP = (pmass - nmass*nmass/pmass)/2.;
  G4ThreeVector const P = nuP*GetRandomVec();
  G4DynamicParticle *nu = 
    new G4DynamicParticle(G4NeutrinoMu::Definition(), P);
  products->PushProducts(nu);
  G4DynamicParticle *neutron = 
    new G4DynamicParticle(G4Neutron::Definition(), -P);
  products->PushProducts(neutron);

  return products;
}


void G4MuPCaptureChannel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,4) != "mu_P" ){
    //    std::ostringstream o;
    std::ostrstream o;
    o<< "Channel " << GetChannelName() << " is not applicable to particle "
     << part->GetParticleName();
    G4Exception("G4MuPCaptureChannel::CheckIsApplicable",
                "MUPC0001", FatalException, 
                o.str());    

  }
}
