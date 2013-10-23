// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, March 19 2010
// ----------------------------------------------------------------

// FIXME ... how can I remove the branching ratio hard coding?

#include "G4MuHe3CaptureChannels.hh"

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
G4MuHe3ProtonChannel::
G4MuHe3ProtonChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_He3->Proton capture", captureRate, verboseLevel),
  psdc(p->GetParticleName(),1.,4,"nu_mu", "proton", "neutron", "neutron")
{
  CheckIsApplicable();
}

G4MuHe3ProtonChannel::
G4MuHe3ProtonChannel(G4MuAtom const* p, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_He3->Proton capture", 0, verboseLevel),
  psdc(p->GetParticleName(),1.,4,"nu_mu", "proton", "neutron", "neutron"){
  CheckIsApplicable();
  G4int const Z = p->GetAtomicNumber();
  G4int const A = p->GetAtomicMass();
  G4int const iSpin = p->GetPDGiSpin();
  G4MuAtomTable const* table = G4MuAtomTable::GetInstance();
  G4MuAtomCaptureRateModel const* model = table->CaptureRateModel(Z,A);
  rate = model->GetCaptureRate(Z,A,iSpin)*10.*CLHEP::perCent; // Measday
}

G4DecayProducts* 
G4MuHe3ProtonChannel::CaptureIt(G4DynamicParticle const* pParticle){
  G4DecayProducts *products =
    psdc.DecayIt(pParticle->GetDefinition()->GetPDGMass());
  products->SetParentParticle(*pParticle);
  return products;
}

void G4MuHe3ProtonChannel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,6) != "mu_He3" ){
#if G4VERSION<4095
    std::ostrstream ed;
#else
    G4ExceptionDescription ed;
#endif
    ed << "Channel " << GetChannelName() << " is not applicable to particle "
       << part->GetParticleName();
    G4Exception("G4MuHe3ProtonChannel::CheckIsApplicable",
                "MUHE30001", FatalException, 
#if G4VERSION<4095
                ed.str());
#else
                ed);
#endif
  }
}

// Deuteron Channel
G4MuHe3DeuteronChannel::
G4MuHe3DeuteronChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_He3->Deuteron capture", captureRate, verboseLevel),
  psdc(p->GetParticleName(),1.,3,"nu_mu", "deuteron", "neutron") {
  CheckIsApplicable();
}

G4MuHe3DeuteronChannel::
G4MuHe3DeuteronChannel(G4MuAtom const* p, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_He3->Deuteron capture", 0, verboseLevel),
  psdc(p->GetParticleName(),1.,3,"nu_mu", "deuteron", "neutron") {
  CheckIsApplicable();
  G4int const Z = p->GetAtomicNumber();
  G4int const A = p->GetAtomicMass();
  G4int const iSpin = p->GetPDGiSpin();
  G4MuAtomTable const* table = G4MuAtomTable::GetInstance();
  G4MuAtomCaptureRateModel const* model = table->CaptureRateModel(Z,A);
  rate = model->GetCaptureRate(Z,A,iSpin)*20.*CLHEP::perCent; // Measday
}

G4DecayProducts* 
G4MuHe3DeuteronChannel::CaptureIt(G4DynamicParticle const* pParticle){
  G4DecayProducts *products =
    psdc.DecayIt(pParticle->GetDefinition()->GetPDGMass());
  products->SetParentParticle(*pParticle);
  return products;
}

void G4MuHe3DeuteronChannel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,6) != "mu_He3" ){
#if G4VERSION<4095
    std::ostrstream ed;
#else
    G4ExceptionDescription ed;
#endif
    ed << "Channel " << GetChannelName() << " is not applicable to particle "
     << part->GetParticleName();
    G4Exception("G4MuHe3DeuteronChannel::CheckIsApplicable",
                "MUHE30002", FatalException, 
#if G4VERSION<4095
                ed.str());
#else
                ed);
#endif
  }
}

// Triton Channel
G4MuHe3TritonChannel::
G4MuHe3TritonChannel(G4MuAtom const* p, G4double captureRate, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_He3->Triton capture", captureRate, verboseLevel),
    psdc(p->GetParticleName(),1.,2,"nu_mu", "triton") {
  CheckIsApplicable();
}

G4MuHe3TritonChannel::
G4MuHe3TritonChannel(G4MuAtom const* p, G4int verboseLevel) :
  G4VMuAtomCaptureKineticsChannel(p, "mu_He3->Triton capture", 0, verboseLevel),
    psdc(p->GetParticleName(),1.,2,"nu_mu", "triton") {
  CheckIsApplicable();
  G4int const Z = p->GetAtomicNumber();
  G4int const A = p->GetAtomicMass();
  G4int const iSpin = p->GetPDGiSpin();
  G4MuAtomTable const* table = G4MuAtomTable::GetInstance();
  G4MuAtomCaptureRateModel const* model = table->CaptureRateModel(Z,A);
  rate = model->GetCaptureRate(Z,A,iSpin)*70.*CLHEP::perCent; // Measday
}

G4DecayProducts* 
G4MuHe3TritonChannel::CaptureIt(G4DynamicParticle const* pParticle){
  G4DecayProducts *products =
    psdc.DecayIt(pParticle->GetDefinition()->GetPDGMass());
  products->SetParentParticle(*pParticle);
  return products;
}

void G4MuHe3TritonChannel::CheckIsApplicable() const {
  if( part->GetParticleName().substr(0,6) != "mu_He3" ){
#if G4VERSION<4095
    std::ostrstream ed;
#else
    G4ExceptionDescription ed;
#endif
    ed << "Channel " << GetChannelName() << " is not applicable to particle "
       << part->GetParticleName();
    G4Exception("G4MuHe3TritonChannel::CheckIsApplicable",
                "MUHE30003", FatalException, 
#if G4VERSION<4095
                ed.str());
#else
                ed);
#endif
  }
}
