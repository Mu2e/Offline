//
// Define a sensitive detector for virtual detectors 
// 
// Original author Ivan Logashenko
//

#include <cstdio>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e incldues
#include "Mu2eG4/inc/CaloReadoutSD.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  CaloReadoutSD::CaloReadoutSD(G4String name, const SimpleConfig& config, CaloCrystalSD *sd)
    : G4VSensitiveDetector(name) {

    G4String HCname("CaloROCollection");
    collectionName.insert(HCname);
    
    crystalSD = sd;

  }


  CaloReadoutSD::~CaloReadoutSD(){ }

  void CaloReadoutSD::Initialize(G4HCofThisEvent* HCE){

    GeomHandle<Calorimeter> cg;
    _nro  = cg->nROPerCrystal();
    _minE = cg->getElectronEmin();

  }
  

  G4bool CaloReadoutSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    //G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();

    const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();

    // Only handle charged events with kinetic energy > 0.1 MeV
    if( aStep->GetTrack()->GetDefinition()->GetPDGCharge() == 0 ) return false;
    if( aStep->GetTrack()->GetKineticEnergy() < _minE ) return false;

    // Get readout ID
    int idro = touchableHandle->GetCopyNumber(0) + touchableHandle->GetCopyNumber(1)*_nro;

    // Ask CaloCrystalSD to create hit

    G4double time = aStep->GetPreStepPoint()->GetGlobalTime();
    G4double edep = aStep->GetTotalEnergyDeposit();

    crystalSD->AddReadoutHit(aStep, idro, time, edep);

    /*
    cout << "CaloReadoutSD: copyNo0=" << touchableHandle->GetCopyNumber(0)
	 << " copyNo1=" << touchableHandle->GetCopyNumber(1)
	 << " time=" << time << " edep=" << edep 
	 << " idro=" << idro 
	 << " charge=" << aStep->GetTrack()->GetDefinition()->GetPDGCharge()
	 << endl;
    */

    return true;

  }

  void CaloReadoutSD::EndOfEvent(G4HCofThisEvent*){
  }

} //namespace mu2e
