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
#include "Mu2eG4/inc/CaloCrystalSD.hh"
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

  G4ThreeVector CaloCrystalSD::_mu2eOrigin;

  CaloCrystalSD::CaloCrystalSD(G4String name, const SimpleConfig& config) :G4VSensitiveDetector(name){
    G4String HCname("CaloCollection");
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


  CaloCrystalSD::~CaloCrystalSD(){ }

  void CaloCrystalSD::Initialize(G4HCofThisEvent* HCE){

    _collection = new StepPointG4Collection
      (SensitiveDetectorName,collectionName[0]); 
    static G4int HCID = -1;
    if(HCID<0){ 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); 
    }
    HCE->AddHitsCollection( HCID, _collection ); 

    _currentSize=0;

  }
  

  G4bool CaloCrystalSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    // Calculate energy deposition in the crystal
    G4double edep = aStep->GetTotalEnergyDeposit();
    if( edep<=0 ) return false;

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
	edm::LogWarning("G4") << "Maximum number of particles reached in CaloCrystalSD: " 
			      << _currentSize << endl;
      }
      return false;
    }

    //G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    //G4int eventId = event->GetEventID();
    //G4int trackId = aStep->GetTrack()->GetTrackID();

    const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();

    // Get crystal ID
    G4int copyNo = touchableHandle->GetCopyNumber(2);

    // Originally the hit position was saved in local crystal frame. 
    // Not it is saved in Mu2e frame, hence the following code is
    // commented out. 
    // Calculate enerdy deposition position along the crystal
    // G4AffineTransform const& toLocal = touchableHandle->GetHistory()->GetTopTransform();
    //G4AffineTransform        toWorld = toLocal.Inverse();
    // G4ThreeVector posWorld = aStep->GetPreStepPoint()->GetPosition();
    // G4ThreeVector posLocal = toLocal.TransformPoint(posWorld);

    StepPointG4* newHit = 
      new StepPointG4(aStep->GetTrack()->GetTrackID(), 
		      copyNo,
                      edep,
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

  void CaloCrystalSD::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      edm::LogWarning("G4") << "Total of " << _currentSize 
			    << " calorimeter hits were generated in the event." 
			    << endl
			    << "Only " << _sizeLimit << " are saved in output collection." 
			    << endl;
    }

    if (verboseLevel>0) { 
      G4int NbHits = _collection->entries();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
             << " hits in the calorimeter: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i]->Print();
    } 

  }

} //namespace mu2e
