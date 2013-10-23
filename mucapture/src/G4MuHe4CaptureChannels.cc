// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, March 19 2010
// ----------------------------------------------------------------

// FIXME ... how can I remove the branching ratio hard coding?

#include "G4MuHe4CaptureChannels.hh"

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4NeutrinoMu.hh"
#include "G4Gamma.hh"
#include "G4ThreeVector.hh"
#include "G4MuAtomTable.hh"

#if G4VERSION<4095
#include <strstream>
#endif

// Proton Channel
G4MuHe4ProtonChannel::
G4MuHe4ProtonChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_He4->Proton capture", captureRate, verboseLevel),
  psdc(p->GetParticleName(),1.,5,"nu_mu") {
  CheckIsApplicable();
  psdc.SetDaughter(1, "proton");
  psdc.SetDaughter(2, "neutron");
  psdc.SetDaughter(3, "neutron");
  psdc.SetDaughter(4, "neutron");
}

G4MuHe4ProtonChannel::
G4MuHe4ProtonChannel(G4MuAtom const* p, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_He4->Proton capture", 0, verboseLevel),
  psdc(p->GetParticleName(),1.,5,"nu_mu") {
  CheckIsApplicable();
  psdc.SetDaughter(1, "proton");
  psdc.SetDaughter(2, "neutron");
  psdc.SetDaughter(3, "neutron");
  psdc.SetDaughter(4, "neutron");
  G4int const Z = p->GetAtomicNumber();
  G4int const A = p->GetAtomicMass();
  G4int const iSpin = p->GetPDGiSpin();
  G4MuAtomTable const* table = G4MuAtomTable::GetInstance();
  G4MuAtomCaptureRateModel const* model = table->CaptureRateModel(Z,A);
  rate = model->GetCaptureRate(Z,A,iSpin)*0.25*CLHEP::perCent; // Measday
}

G4DecayProducts* 
G4MuHe4ProtonChannel::CaptureIt(G4DynamicParticle const* pParticle){
  G4DecayProducts *products =
    psdc.DecayIt(pParticle->GetDefinition()->GetPDGMass());
  products->SetParentParticle(*pParticle);
  return products;
}

void G4MuHe4ProtonChannel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,6) != "mu_He4" ){
#if G4VERSION<4095
    std::ostrstream ed;
#else
    G4ExceptionDescription ed;
#endif
    ed << "Channel " << GetChannelName() << " is not applicable to particle "
     << part->GetParticleName();
    G4Exception("G4MuHe4ProtonChannel::CheckIsApplicable",
                "MUHE40001", FatalException, 
#if G4VERSION<4095
                ed.str());
#else
                ed);
#endif
  }
}

// Deuteron Channel
G4MuHe4DeuteronChannel::
G4MuHe4DeuteronChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_He4->Deuteron capture", captureRate, verboseLevel),
    psdc(p->GetParticleName(),1.,4,"nu_mu", "deuteron", "neutron", "neutron") {
  CheckIsApplicable();
} 

G4MuHe4DeuteronChannel::
G4MuHe4DeuteronChannel(G4MuAtom const* p, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_He4->Deuteron capture", 0, verboseLevel),
    psdc(p->GetParticleName(),1.,4,"nu_mu", "deuteron", "neutron", "neutron") {
  CheckIsApplicable();
  G4int const Z = p->GetAtomicNumber();
  G4int const A = p->GetAtomicMass();
  G4int const iSpin = p->GetPDGiSpin();
  G4MuAtomTable const* table = G4MuAtomTable::GetInstance();
  G4MuAtomCaptureRateModel const* model = table->CaptureRateModel(Z,A);
  rate = model->GetCaptureRate(Z,A,iSpin)*2.*CLHEP::perCent; // Measday
}

G4DecayProducts* 
G4MuHe4DeuteronChannel::CaptureIt(G4DynamicParticle const* pParticle){
  G4DecayProducts *products =
    psdc.DecayIt(pParticle->GetDefinition()->GetPDGMass());
  products->SetParentParticle(*pParticle);
  return products;
}

void G4MuHe4DeuteronChannel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,6) != "mu_He4" ){
#if G4VERSION<4095
    std::ostrstream ed;
#else
    G4ExceptionDescription ed;
#endif
    ed << "Channel " << GetChannelName() << " is not applicable to particle "
     << part->GetParticleName();
    G4Exception("G4MuHe4DeuteronChannel::CheckIsApplicable",
                "MUHE40002", FatalException, 
#if G4VERSION<4095
                ed.str());
#else
                ed);
#endif
  }
}

// Triton Channel
G4MuHe4TritonChannel::
G4MuHe4TritonChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_He4->Triton capture", captureRate, verboseLevel),
    psdc(p->GetParticleName(),1.,3,"nu_mu", "triton", "neutron") {
  CheckIsApplicable();
}

G4MuHe4TritonChannel::
G4MuHe4TritonChannel(G4MuAtom const* p, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_He4->Triton capture", 0, verboseLevel),
  psdc(p->GetParticleName(),1.,3,"nu_mu", "triton", "neutron") {
  CheckIsApplicable();
  G4int const Z = p->GetAtomicNumber();
  G4int const A = p->GetAtomicMass();
  G4int const iSpin = p->GetPDGiSpin();
  G4MuAtomTable const* table = G4MuAtomTable::GetInstance();
  G4MuAtomCaptureRateModel const* model = table->CaptureRateModel(Z,A);
  rate = model->GetCaptureRate(Z,A,iSpin)*97.75*CLHEP::perCent; // Measday
}

G4DecayProducts* 
G4MuHe4TritonChannel::CaptureIt(G4DynamicParticle const* pParticle){
  G4DecayProducts *products =
    psdc.DecayIt(pParticle->GetDefinition()->GetPDGMass());
  products->SetParentParticle(*pParticle);
  return products;
}

void G4MuHe4TritonChannel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,6) != "mu_He4" ){
 #if G4VERSION<4095
    std::ostrstream ed;
#else
    G4ExceptionDescription ed;
#endif
    ed<< "Channel " << GetChannelName() << " is not applicable to particle "
     << part->GetParticleName();
    G4Exception("G4MuHe4TritonChannel::CheckIsApplicable",
                "MUHE40003", FatalException, 
#if G4VERSION<4095
                ed.str());
#else
                ed);
#endif
  }
}
