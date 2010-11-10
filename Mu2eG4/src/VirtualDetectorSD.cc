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
#include "Mu2eG4/inc/VirtualDetectorSD.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "Mu2eG4/inc/StepPointG4.hh"

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

  G4ThreeVector VirtualDetectorSD::_mu2eOrigin;

  VirtualDetectorSD::VirtualDetectorSD(G4String name, const SimpleConfig& config) :G4VSensitiveDetector(name){
    G4String HCname("VDCollection");
    collectionName.insert(HCname);

    // Get list of events for which to make debug printout.
    string key("g4.virtualSDEventList");
    if ( config.hasName(key) ){
      vector<int> list;
      config.getVectorInt(key,list);
      _debugList.add(list);
    }

    _sizeLimit = config.getInt("g4.stepsSizeLimit",0);
    _currentSize = 0;

  }


  VirtualDetectorSD::~VirtualDetectorSD(){ }

  void VirtualDetectorSD::Initialize(G4HCofThisEvent* HCE){

    _collection = new StepPointG4Collection
      (SensitiveDetectorName,collectionName[0]); 
    static G4int HCID = -1;
    if(HCID<0){ 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); 
    }
    HCE->AddHitsCollection( HCID, _collection ); 

    _currentSize=0;

  }
  

  G4bool VirtualDetectorSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
	edm::LogWarning("G4") << "Maximum number of particles reached in VirtualDetectorSD: " 
			      << _currentSize << endl;
      }
      return false;
    }

    //G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    //const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();
    //G4int eventId = event->GetEventID();
    //G4int trackId = aStep->GetTrack()->GetTrackID();
    //G4int copyNo = touchableHandle->GetCopyNumber();

    // The points coordinates are saved in the mu2e world

    StepPointG4* newHit = 
      new StepPointG4(aStep->GetTrack()->GetTrackID(),
                      aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetCopyNo(),
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

  void VirtualDetectorSD::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      edm::LogWarning("G4") << "Total of " << _currentSize 
			    << " virtual detector hits were generated in the event." 
			    << endl
			    << "Only " << _sizeLimit << " are saved in output collection." 
			    << endl;
      cout << "Total of " << _currentSize 
	   << " virtual detector hits were generated in the event." 
	   << endl
	   << "Only " << _sizeLimit << " are saved in output collection." 
	   << endl;
    }

    if (verboseLevel>0) { 
      G4int NbHits = _collection->entries();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
             << " hits in the virtual detectors: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i]->Print();
    } 

  }

} //namespace mu2e
