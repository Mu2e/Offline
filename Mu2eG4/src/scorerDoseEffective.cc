#include "Offline/Mu2eG4/inc/scorerDoseEffective.hh"
#include "Offline/Mu2eG4/inc/scorerFTDConverter.hh"

#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"

#include "G4Box.hh"


namespace mu2e {


  scorerDoseEffective::scorerDoseEffective(const G4String& name,G4int depth)
    : G4VPrimitiveScorer(name,depth),
      fDepthi(2),
      fDepthj(1),
      fDepthk(0),
      FTDConverter("ISO")
  {}


  void scorerDoseEffective::Initialize(G4HCofThisEvent* HCE)
  {
     EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
     if(HCID < 0) HCID = GetCollectionID(0);
     HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
  }

  void scorerDoseEffective::clear() { EvtMap->clear(); }


  G4bool scorerDoseEffective::ProcessHits(G4Step* aStep, G4TouchableHistory*)
  {
    G4double stepLength = aStep->GetStepLength();
    if (stepLength <1e-12) return false;

    G4int idx = ((G4TouchableHistory*) (aStep->GetPreStepPoint()->GetTouchable()))
                ->GetReplicaNumber(indexDepth);
    G4double cubicVolume = ComputeVolume(aStep, idx);

    G4double CellFlux = stepLength / cubicVolume;
    if (weighted) CellFlux *= aStep->GetPreStepPoint()->GetWeight();


    //table lookup
    int pdgCode      = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    double energy    = aStep->GetPreStepPoint()->GetKineticEnergy();
    double FTDfactor = FTDConverter.evaluate(pdgCode,energy);
    CellFlux *= FTDfactor;

    G4int index = GetIndex(aStep);
    EvtMap->add(index, CellFlux);

    return true;
  }


  /*
  G4bool scorerDose::ProcessHits(G4Step* aStep,G4TouchableHistory*)
  {
    G4StepPoint* preStep = aStep->GetPreStepPoint();
    G4StepPoint* posStep = aStep->GetPostStepPoint();
    G4bool IsEnter = preStep->GetStepStatus()==fGeomBoundary;
    G4bool IsExit  = posStep->GetStepStatus()==fGeomBoundary;

    G4int index = GetIndex(aStep);

    G4bool flag = FALSE;
    if ( IsEnter && fDirection == fCurrent_In ) flag = TRUE;
    else if ( IsExit  && fDirection == fCurrent_Out ) flag = TRUE;
    else if ( (IsExit||IsEnter) && fDirection == fCurrent_InOut  ) flag = TRUE;

    if ( flag ){
      G4double val = 1.0;
      if (weighted) val *= aStep->GetPreStepPoint()->GetWeight();
      //EvtMap->add(index,val);
    }

     //this option takes the actual physica volume - ok if volumes are smaller than mesh size
    G4VPhysicalVolume* physVol = aStep->GetTrack()->GetStep()->GetPreStepPoint()->GetPhysicalVolume();
    //this option takes the mesh size (same as DoseDepoisted scorer), ok if volumes are much larger than mesh size
    //G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
    G4VPVParameterisation* physParam = physVol->GetParameterisation();
    G4VSolid* solid = 0;
    if (physParam) {
      if (index>0) {
        solid = physParam->ComputeSolid(index, physVol);
        solid->ComputeDimensions(physParam,index,physVol);
      }
    }
    else solid = physVol->GetLogicalVolume()->GetSolid();

    double cubicVolume = (solid) ?  solid->GetCubicVolume() : 0;
    G4double edep = aStep->GetTotalEnergyDeposit();
    G4double density = aStep->GetTrack()->GetStep()->GetPreStepPoint()->GetMaterial()->GetDensity();
    G4double dose    = cubicVolume>1e-20 ? edep / ( density * cubicVolume )/CLHEP::gray : 0;  // unit 1e-12 Gy

    //dumpStep(aStep,solid);

    if ( edep <1e-20 ) return FALSE;
    EvtMap->add(index,dose);
    return TRUE;
  }
  */

  G4double scorerDoseEffective::ComputeVolume(G4Step* aStep, G4int idx)
  {
     G4VSolid* solid = ComputeSolid(aStep, idx);
     assert(solid);
     return solid->GetCubicVolume();
  }

  G4int scorerDoseEffective::GetIndex(G4Step* aStep)
  {
    const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int i = touchable->GetReplicaNumber(fDepthi);
    G4int j = touchable->GetReplicaNumber(fDepthj);
    G4int k = touchable->GetReplicaNumber(fDepthk);
    return i*fNj*fNk+j*fNk+k;
  }


  void scorerDoseEffective::PrintAll()
  {
    G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
    G4cout << " PrimitiveScorer " << GetName() << G4endl;
    G4cout << " Number of entries " << EvtMap->entries() << G4endl;
    std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
    for(; itr != EvtMap->GetMap()->end(); itr++) {
      G4cout << "  copy no.: " << itr->first
             << "  track count: " << *(itr->second)
             << " [tracks] "
             << G4endl;
    }
  }


}

