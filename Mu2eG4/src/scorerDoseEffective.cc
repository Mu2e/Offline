#include "Offline/Mu2eG4/inc/scorerDoseEffective.hh"
#include "Offline/Mu2eG4/inc/scorerFTDConverter.hh"
#include "CLHEP/Random/RandFlat.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"

#include "G4Box.hh"



namespace mu2e {
  scorerDoseEffective::scorerDoseEffective(const G4String& name,
                                           const Mu2eG4Config::Physics& configPhysics,
                                           G4int depth) :
    G4VPrimitiveScorer(name,depth),
    fDepthi_(2),
    fDepthj_(1),
    fDepthk_(0),
    FTDConverter_(configPhysics.radiationTableName())
  {}


  void scorerDoseEffective::Initialize(G4HCofThisEvent* HCE)
  {
    EvtMap_ = new G4THitsMap<G4double>(detector->GetName(), GetName());
    if (HCID_ < 0) HCID_ = GetCollectionID(0);
    HCE->AddHitsCollection(HCID_, (G4VHitsCollection*)EvtMap_);
  }

  void scorerDoseEffective::clear() { EvtMap_->clear(); }


  G4bool scorerDoseEffective::ProcessHits(G4Step* aStep, G4TouchableHistory*)
  {
    G4double stepLength = aStep->GetStepLength();
    //if (stepLength <1e-12) return false;

    G4int idx  = ((G4TouchableHistory*) (aStep->GetPreStepPoint()->GetTouchable()))
                  ->GetReplicaNumber(indexDepth);
    G4double cubicVolume = ComputeVolume(aStep, idx);
    G4double CellFlux    = stepLength / cubicVolume;
    CellFlux            *= aStep->GetPreStepPoint()->GetWeight();


    //table lookup
    int pdgCode      = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    double energy    = aStep->GetPreStepPoint()->GetKineticEnergy();
    double FTDfactor = FTDConverter_.evaluate(pdgCode,energy);
    CellFlux *= FTDfactor;

    G4int index = GetIndex(aStep);
    EvtMap_->add(index, CellFlux);

    return true;
  }

  G4double scorerDoseEffective::ComputeVolume(G4Step* aStep, G4int idx)
  {
    G4VSolid* solid = ComputeSolid(aStep, idx);
    assert(solid);
    return solid->GetCubicVolume();
  }

  G4int scorerDoseEffective::GetIndex(G4Step* aStep)
  {
    const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int i = touchable->GetReplicaNumber(fDepthi_);
    G4int j = touchable->GetReplicaNumber(fDepthj_);
    G4int k = touchable->GetReplicaNumber(fDepthk_);
    return i*fNj*fNk+j*fNk+k;
  }


  void scorerDoseEffective::PrintAll()
  {
    G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
    G4cout << " PrimitiveScorer " << GetName() << G4endl;
    G4cout << " Number of entries " << EvtMap_->entries() << G4endl;
    std::map<G4int,G4double*>::iterator itr = EvtMap_->GetMap()->begin();
    for(; itr != EvtMap_->GetMap()->end(); itr++) {
      G4cout << "  copy no.: " << itr->first
             << "  track count: " << *(itr->second)
             << " [tracks] "
             << G4endl;
    }
  }


}

