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

  G4ThreeVector CaloReadoutSD::_mu2eOrigin;

  CaloReadoutSD::CaloReadoutSD(G4String name, const SimpleConfig& config)
    : G4VSensitiveDetector(name) {

    G4String HCname("CaloROCollection");
    collectionName.insert(HCname);
    
    // Get list of events for which to make debug printout.
    string key("g4.calorimeterSDEventList");
    if ( config.hasName(key) ){
      vector<int> list;
      config.getVectorInt(key,list);
      _debugList.add(list);
    }

    _sizeLimit = config.getInt("g4.stepsSizeLimit",0);

  }


  CaloReadoutSD::~CaloReadoutSD(){ }

  void CaloReadoutSD::Initialize(G4HCofThisEvent* HCE){

    _collection = new StepPointG4Collection
      (SensitiveDetectorName,collectionName[0]); 
    static G4int HCID = -1;
    if(HCID<0){ 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); 
    }
    HCE->AddHitsCollection( HCID, _collection ); 

    _currentSize=0;

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

    // Check that number of steps did not exceed the limit 

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
	edm::LogWarning("G4") << "Maximum number of particles reached in CaloCrystalSD: " 
			      << _currentSize << endl;
      }
      return false;
    }

    // Get readout ID
    int idro = touchableHandle->GetCopyNumber(0) + touchableHandle->GetCopyNumber(1)*_nro;

    // The points coordinates are saved in the mu2e world

    StepPointG4* newHit = 
      new StepPointG4(aStep->GetTrack()->GetTrackID(), 
		      idro,
                      aStep->GetTotalEnergyDeposit(),
                      aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin,
                      aStep->GetPreStepPoint()->GetMomentum(),
                      aStep->GetPreStepPoint()->GetGlobalTime(),
                      aStep->GetPreStepPoint()->GetProperTime(),
                      aStep->GetStepLength()
                      );

    // The collection takes ownership of the hit. 
    _collection->insert( newHit );

    return true;

  }

  void CaloReadoutSD::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      edm::LogWarning("G4") << "Total of " << _currentSize 
			    << " calorimeter RO hits were generated in the event." 
			    << endl
			    << "Only " << _sizeLimit << " are saved in output collection." 
			    << endl;
    }

    if (verboseLevel>0) { 
      G4int NbHits = _collection->entries();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
             << " RO hits in the calorimeter: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i]->Print();
    } 

  }

} //namespace mu2e
