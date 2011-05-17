//
// Defines sensitive detector for CRSScintillatorBar
//
// $Id: CRSScintillatorBarSD.cc,v 1.2 2011/05/17 15:36:00 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:36:00 $
//
// Original author KLG 
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e incldues
#include "Mu2eG4/inc/CRSScintillatorBarSD.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  G4ThreeVector CRSScintillatorBarSD::_mu2eOrigin;

  CRSScintillatorBarSD::CRSScintillatorBarSD(G4String name, const SimpleConfig& config):
    G4VSensitiveDetector(name),
    _collection(0),
    _debugList(0),
    _sizeLimit(config.get<int>("g4.stepsSizeLimit",0)),
    _currentSize(0)
  {

    // Get list of events for which to make debug printout.
    string key("g4.CRSScintillatorBarSDEventList");
    if ( config.hasName(key) ){
      vector<int> list;
      config.getVectorInt(key,list);
      _debugList.add(list);
    }

  }


  CRSScintillatorBarSD::~CRSScintillatorBarSD(){ }

  void CRSScintillatorBarSD::Initialize(G4HCofThisEvent* HCE){

    _currentSize=0;

  }
  

  G4bool CRSScintillatorBarSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
	mf::LogWarning("G4") << "Maximum number of particles reached in CRSScintillatorBarSD: " 
			      << _currentSize << endl;
      }
      return false;
    }

    // The points coordinates are saved relative to mu2eOrigin 

    // we add the hit to the framework collection
    _collection->
      push_back(StepPointMC(aStep->GetTrack()->GetTrackID(),
                            aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetCopyNo(),
                            aStep->GetTotalEnergyDeposit(),
                            aStep->GetPreStepPoint()->GetGlobalTime(),
                            aStep->GetPreStepPoint()->GetProperTime(),
                            aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin,
                            aStep->GetPreStepPoint()->GetMomentum(),
                            aStep->GetStepLength()
                            ));

    return true;

  }


  void CRSScintillatorBarSD::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      mf::LogWarning("G4") << "Total of " << _currentSize 
			    << " CRS Scintillator Bar hits were generated in the event." 
			    << endl
			    << "Only " << _sizeLimit << " are saved in output collection." 
			    << endl;
      cout << "Total of " << _currentSize 
	   << " CRS Scintillator Bar hits were generated in the event." 
	   << endl
	   << "Only " << _sizeLimit << " are saved in output collection." 
	   << endl;
    }

    if (verboseLevel>0) { 
      G4int NbHits = _collection->size();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
             << " hits in the CRS Scintillator Bar: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i].print(G4cout);
    } 

  }


  void CRSScintillatorBarSD::beforeG4Event(StepPointMCCollection& outputHits) {
    _collection = &outputHits;
    return;
  } // end of beforeG4Event


} //namespace mu2e
