// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

// FIXME ... need lots of verbose output code

#include "G4MuonMinusAtomicCapture.hh"
#include "G4MuonMinus.hh"
#include "G4MuAtom.hh"
#include "G4MuAtomTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4HadronicProcessType.hh"
#include "G4HadronicProcessStore.hh"

// CLHEP includes
#include "CLHEP/Units/PhysicalConstants.h"

G4MuonMinusAtomicCapture::  
G4MuonMinusAtomicCapture(G4String const& name, G4ProcessType type) : 
  G4VRestProcess(name,type),
  pSelector(new G4ElementSelector)
{
  SetProcessSubType(fCapture); // FIXME ... possibly need to register a new process
                               // subtype with someone...
}


G4MuonMinusAtomicCapture::~G4MuonMinusAtomicCapture(){
  delete pSelector;
}

G4bool G4MuonMinusAtomicCapture::IsApplicable(G4ParticleDefinition const& p){
  return ( &p == G4MuonMinus::Definition() );
}

G4VParticleChange* 
G4MuonMinusAtomicCapture::AtRestDoIt(const G4Track& track, const G4Step& step){ 
  aParticleChange.Initialize(track);

  // select element and get Z,A.
  // G4Element const* aEle = pSelector->GetElement(track.GetMaterial());

  // replacing deprecated G4StopElementSelector with G4ElementSelector
  // the fuctionality below is probably a duplication of what is in
  // the G4ElementSelector

  G4Nucleus* nucleus = new G4Nucleus(); // the next call sets its Z & A, we will ignore it
  G4Element const* aEle = pSelector->SelectZandA(track, nucleus);

  G4int const targetZ = aEle->GetZ();
  G4int targetA;
  G4int ni = 0;

  G4IsotopeVector const* isv = aEle->GetIsotopeVector();
  if(isv) 
    ni = isv->size();

  if(ni == 0){
    targetA = aEle->GetN();
  } else if(ni == 1) {
    targetA = aEle->GetIsotope(0)->GetN();
  } else {
    G4double const* ab = aEle->GetRelativeAbundanceVector();
    G4double y = G4UniformRand();
    G4int j = -1;
    ni--;
    do {
      j++;
      y -= ab[j];
    } while (y > 0.0 && j < ni);
    targetA = aEle->GetIsotope(j)->GetN();
  }

  verboseLevel >=1 &&
    G4cout << "targetZ: " << targetZ << " targetA: " << targetA << '\n';

  G4int nSecondaries = 0;

  G4double const globalTime = track.GetGlobalTime();
  G4ThreeVector const position = track.GetPosition();

  // Now, we've got the target ... we need to build the MuAtom 
  G4DynamicParticle *dynMuAtom = new G4DynamicParticle();
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  G4int iSpin = table->StateModel(targetZ, targetA)->GetSpinState(step);
  G4MuAtom *muAtom = table->GetMuAtom(targetZ, targetA, iSpin); // FIXME this leaks
  dynMuAtom->SetDefinition(muAtom);
  dynMuAtom->SetMomentum(G4ThreeVector()); // how is the muon energy dissipated ???
  G4int charge = table->ChargeModel(targetZ, targetA)->GetCharge(step); // this is also channel insertion
  dynMuAtom->SetCharge(charge); // sets in units of eplus
  // FIXME ... need to include the electron binding energy contribution...
  dynMuAtom->SetMass(CLHEP::electron_mass_c2*(targetZ-1)+muAtom->GetPDGMass());
  G4Track* muAtomTrack = new G4Track(dynMuAtom, globalTime, position);
  muAtomTrack->SetTouchableHandle( track.GetTouchableHandle() );
  ++nSecondaries;

  // Run the electromagnetic cascade model
  
  // if( G4MuonMinusAtomicCaptureCascadeModel * model = 
  //     table->CascadeModel(targetZ, targetA) ){
  //   // IMPLEMENT_ME ???
  // }

  // it looks like the entire em cascade, capture, decay of bound muon is missing???
  // but what is in G4MuAtomGenericCaptureChannel???


  aParticleChange.SetNumberOfSecondaries(nSecondaries);
  aParticleChange.AddSecondary( muAtomTrack );

  // the above "converts" mu+atom into muatom...

  // Add the EM secondaries from the cascade 

  aParticleChange.ProposeLocalEnergyDeposit(0.0);
  aParticleChange.ProposeTrackStatus(fStopAndKill); // we kill the muon and create the muAtom

  return &aParticleChange;
}


// FIXME do we need those???
void G4MuonMinusAtomicCapture::PreparePhysicsTable(const G4ParticleDefinition& p) 
{
  G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this, &p);
}

void G4MuonMinusAtomicCapture::BuildPhysicsTable(const G4ParticleDefinition& p) 
{
  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}



