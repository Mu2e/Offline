#include "Offline/Mu2eG4/inc/customScorer.hh"
#include "G4UnitsTable.hh"
#include "G4VPrimitiveScorer.hh"

#include "G4VSolid.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "G4UnitsTable.hh"

namespace mu2e {

  customScorer::customScorer(G4String name, G4int direction, G4int depth)
    : G4VPrimitiveScorer(name,depth),HCID(-1),fDirection(direction),weighted(false), fDepthi(2),fDepthj(1),fDepthk(0)
  {
    fNi=1;
    fNj=1;
    fNk=1;
    SetUnit("");
  }

  G4bool customScorer::ProcessHits(G4Step* aStep,G4TouchableHistory*)
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



  void customScorer::Initialize(G4HCofThisEvent* HCE)
  {
    EvtMap = new G4THitsMap<G4double>(detector->GetName(),GetName());
    if(HCID < 0) {HCID = GetCollectionID(0);}
    HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
    G4cout <<   "----- PS Initialize " << G4endl;
  }


  void customScorer::clear(){
    EvtMap->clear();
  }

  void customScorer::PrintAll()
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

  void customScorer::SetUnit(const G4String& unit)
  {
    if (unit == "" ){
      unitName = unit;
      unitValue = 1.0;
    } else {
        G4String msg = "Invalid unit ["+unit+"] (Current  unit is [" +GetUnit()+"] ) for " + GetName();
        G4Exception("customScorer::SetUnit","DetPS0018",JustWarning,msg);
    }

  }

  G4int customScorer::GetIndex(G4Step* aStep)
  {
    const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int i = touchable->GetReplicaNumber(fDepthi);
    G4int j = touchable->GetReplicaNumber(fDepthj);
    G4int k = touchable->GetReplicaNumber(fDepthk);
    return i*fNj*fNk+j*fNk+k;
  }


  void customScorer::dumpStep(G4Step* aStep, G4VSolid* solid)
  {
    const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();
    G4ThreeVector posWorld = aStep->GetPreStepPoint()->GetPosition();
    for (int i=0;i<=touchableHandle->GetHistoryDepth();++i)
      std::cout<<"Transform level "<<i<<"   "<<touchableHandle->GetCopyNumber(i)
               <<"   "<<touchableHandle->GetHistory()->GetTransform(touchableHandle->GetHistoryDepth()-i).TransformPoint(posWorld)
               <<"  "<<touchableHandle->GetVolume(i)->GetTranslation()
               <<"  "<<touchableHandle->GetVolume(i)->GetName()<<std::endl;

    G4cout << " trk " << aStep->GetTrack()->GetTrackID()
           << " Momentum " <<aStep->GetPreStepPoint()->GetMomentum()
           << " ID " << aStep->GetTrack()->GetParticleDefinition()->GetParticleName()
           <<    G4endl;

    double cubicVolume = (solid)?  solid->GetCubicVolume() : 0;
    G4double edep      = aStep->GetTotalEnergyDeposit();
    G4double density   = aStep->GetTrack()->GetStep()->GetPreStepPoint()->GetMaterial()->GetDensity();
    G4double dose      = edep/ ( density * cubicVolume );  // unit 1e-12 Gy

    G4cout<<"Solid "<<solid->GetName()<<" Vol="<<cubicVolume<<"  Density="<<density<<" Dose="<<dose<<std::endl;
    G4cout<<dose/CLHEP::gray<<std::endl;

    G4cout<<std::endl;
  }

}
