// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, March 18 2010
// ----------------------------------------------------------------

#include "G4MuDCaptureChannel.hh"

#include "G4Neutron.hh"
#include "G4NeutrinoMu.hh"
#include "G4ThreeVector.hh"
#include "G4MuAtomTable.hh"

#include <sstream>
#include <strstream>

// base class

G4MuDCaptureChannel::
G4MuDCaptureChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_D capture", captureRate, verboseLevel),
  psdc(p->GetParticleName(),1.,3,"nu_mu","neutron","neutron") {
  CheckIsApplicable();
}

G4MuDCaptureChannel::
G4MuDCaptureChannel(G4MuAtom const* p, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_D capture", 0., verboseLevel),
  psdc(p->GetParticleName(),1.,3,"nu_mu","neutron","neutron") {
  CheckIsApplicable();
  G4int const Z = p->GetAtomicNumber();
  G4int const A = p->GetAtomicMass();
  G4int const iSpin = p->GetPDGiSpin();
  G4MuAtomTable const* table = G4MuAtomTable::GetInstance();
  G4MuAtomCaptureRateModel const* model = table->CaptureRateModel(Z,A);
  rate = model->GetCaptureRate(Z,A,iSpin);
}

G4DecayProducts* G4MuDCaptureChannel::CaptureIt(G4DynamicParticle const* pParticle){
  // leverage existing code...
  G4DecayProducts *products = 
    psdc.DecayIt(pParticle->GetDefinition()->GetPDGMass());
  products->SetParentParticle(*pParticle);
  return products;
}

void G4MuDCaptureChannel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,4) != "mu_D" ){
    //    std::ostringstream o;
    std::ostrstream o;
    o<< "Channel " << GetChannelName() << " is not applicable to particle "
     << part->GetParticleName();
    G4Exception("G4MuDCaptureChannel::CheckIsApplicable",
                "MUDC0001", FatalException, 
                o.str());
  }
}
