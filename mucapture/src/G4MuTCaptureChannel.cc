// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, March 18 2010
// ----------------------------------------------------------------

#include "G4MuTCaptureChannel.hh"

#include "G4Neutron.hh"
#include "G4NeutrinoMu.hh"
#include "G4ThreeVector.hh"
#include "G4MuAtomTable.hh"

#include <strstream>

// base class

G4MuTCaptureChannel::
G4MuTCaptureChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_T capture", captureRate, verboseLevel),
    psdc(p->GetParticleName(),1.,4,"nu_mu","neutron","neutron","neutron") {
  CheckIsApplicable();
}

G4MuTCaptureChannel::
G4MuTCaptureChannel(G4MuAtom const* p, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_T capture", 0, verboseLevel),
    psdc(p->GetParticleName(),1.,4,"nu_mu","neutron","neutron","neutron") {
  CheckIsApplicable();
  G4int const Z = p->GetAtomicNumber();
  G4int const A = p->GetAtomicMass();
  G4int const iSpin = p->GetPDGiSpin();
  G4MuAtomTable const* table = G4MuAtomTable::GetInstance();
  G4MuAtomCaptureRateModel const* model = table->CaptureRateModel(Z,A);
  rate = model->GetCaptureRate(Z,A,iSpin);
}

G4DecayProducts* G4MuTCaptureChannel::CaptureIt(G4DynamicParticle const* pParticle){
  G4DecayProducts *products = 
    psdc.DecayIt(pParticle->GetDefinition()->GetPDGMass());
  products->SetParentParticle(*pParticle);
  return products;
}

void G4MuTCaptureChannel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,4) != "mu_T" ){
    //    std::ostringstream o;
    std::ostrstream o;
    o<< "Channel " << GetChannelName() << " is not applicable to particle "
     << part->GetParticleName();
    G4Exception("G4MuTCaptureChannel::CheckIsApplicable",
                "MUPT0001", FatalException, 
                o.str());    
  }
}
