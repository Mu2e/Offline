#include "Offline/Mu2eG4/inc/scorerDelayedDose.hh"
#include "Offline/Mu2eG4/inc/scorerFTDConverter.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"



namespace mu2e {
  scorerDelayedDose::scorerDelayedDose(const G4String& name,
                                       const Mu2eG4Config::Physics& configPhysics,
                                       G4int depth) :
    G4VPrimitiveScorer(name,depth),
    fDepthi_       (2),
    fDepthj_       (1),
    fDepthk_       (0),
    FTDConverter_  (configPhysics.radiationTableName())
  {
    readTimeProfile(configPhysics.coolTimeProfileRad());
  }

  void scorerDelayedDose::readTimeProfile(std::string filename)
  {
    bin_.clear();
    profile_.clear();

    std::ifstream infile (filename, std::ios::in);
    if (!infile){
       throw cet::exception("INIT")<<"scorerDelayedDose::readTimeProfile "
                                   <<filename<<" does not exist\n";
    }

    G4double bin(0.), flux(0.),rsum(0.);
    while (infile >> bin >> flux) {
      bin_.push_back(bin * s);
      rsum += flux;
      profile_.push_back(flux);
    }

    if (bin_.size()<2) {
       throw cet::exception("INIT")<<"scorerDelayedDose::readTimeProfile "
                                   <<filename<<" has wrong format, less than 2 entries\n";
    }

    for (size_t i=1;i<bin_.size();++i){
      if (profile_[i-1]>1e-6) totalCoolTime_ += (bin_[i]-bin_[i-1])/s;
    }
  }

  void scorerDelayedDose::Initialize(G4HCofThisEvent* HCE)
  {
    EvtMap_ = new G4THitsMap<G4double>(detector->GetName(), GetName());
    if (HCID_ < 0) HCID_ = GetCollectionID(0);
    HCE->AddHitsCollection(HCID_, (G4VHitsCollection*)EvtMap_);
  }

  void scorerDelayedDose::clear() { EvtMap_->clear(); }


  G4bool scorerDelayedDose::ProcessHits(G4Step* aStep, G4TouchableHistory*)
  {
    G4double stepLength = aStep->GetStepLength();
    //if (stepLength <1e-12) return false;

    if (!IsInDecayWindow(aStep->GetPreStepPoint()->GetGlobalTime())) return false;

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

    //Normalize the flux in seconds
    CellFlux /= totalCoolTime_;

    G4int index = GetIndex(aStep);
    EvtMap_->add(index, CellFlux);

    return true;
  }


  G4bool scorerDelayedDose::IsInDecayWindow(G4double time)
  {
    for (size_t i=1;i<bin_.size();++i){
      if (time > bin_[i-1] && time < bin_[i] && profile_[i-1]>1e-6) return true;
    }
    return false;
  }

  G4double scorerDelayedDose::ComputeVolume(G4Step* aStep, G4int idx)
  {
    G4VSolid* solid = ComputeSolid(aStep, idx);
    assert(solid);
    return solid->GetCubicVolume();
  }

  G4int scorerDelayedDose::GetIndex(G4Step* aStep)
  {
    const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int i = touchable->GetReplicaNumber(fDepthi_);
    G4int j = touchable->GetReplicaNumber(fDepthj_);
    G4int k = touchable->GetReplicaNumber(fDepthk_);
    return i*fNj*fNk+j*fNk+k;
  }


  void scorerDelayedDose::PrintAll()
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
