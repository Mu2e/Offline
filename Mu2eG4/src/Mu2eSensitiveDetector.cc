//
// Defines sensitive detector for a typicaly numbered volume using mu2e reference frame
//
// $Id: Mu2eSensitiveDetector.cc,v 1.1 2012/05/29 22:53:01 genser Exp $
// $Author: genser $
// $Date: 2012/05/29 22:53:01 $
//
// Original author KLG
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  Mu2eSensitiveDetector::Mu2eSensitiveDetector(G4String const name, const SimpleConfig& config):
    G4VSensitiveDetector(name),
    _collection(0),
    _processInfo(0),
    _mu2eOrigin(GeomHandle<WorldG4>()->mu2eOriginInWorld()),
    _debugList(0),
    _sizeLimit(config.getInt("g4.stepsSizeLimit",0)),
    _currentSize(0),
    _simID(0),
    _event(0)
  {

   // Get list of events for which to make debug printout.
    // we generate the key string based on the detector name
    // consult Mu2eG4/inc/SensitiveDetectorName.hh for the names

    std::ostringstream sdKeyName;
    sdKeyName<<"g4."<< SensitiveDetectorName << "SDEventList";
    //G4cout << __func__ << " sdKeyName: " << sdKeyName.str() << G4endl;
 
    string key(sdKeyName.str());
    if ( config.hasName(key) ){
      vector<int> list;
      config.getVectorInt(key,list);
      _debugList.add(list);
    }

  }

  void Mu2eSensitiveDetector::Initialize(G4HCofThisEvent* HCE){

    _currentSize=0;

  }


  G4bool Mu2eSensitiveDetector::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    _currentSize += 1;

    if ( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in " 
                             << SensitiveDetectorName
                             << ": "
                             << _currentSize << endl;
      }
      return false;
    }

    // Which process caused this step to end?
    G4String const& pname  = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    ProcessCode endCode(_processInfo->findAndCount(pname));

    // Add the hit to the framework collection.
    // The point's coordinates are saved in the mu2e coordinate system.
    _collection->
      push_back(StepPointMC(art::Ptr<SimParticle>
                            ( *_simID, 
                              aStep->GetTrack()->GetTrackID(), 
                              _event->productGetter(*_simID) ),
                            aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(),
                            aStep->GetTotalEnergyDeposit(),
                            aStep->GetNonIonizingEnergyDeposit(),
                            aStep->GetPreStepPoint()->GetGlobalTime(),
                            aStep->GetPreStepPoint()->GetProperTime(),
                            aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin,
                            aStep->GetPreStepPoint()->GetMomentum(),
                            aStep->GetStepLength(),
                            endCode
                            ));

    return true;

  }


  void Mu2eSensitiveDetector::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      mf::LogWarning("G4") << "Total of " << _currentSize << " "
                           << SensitiveDetectorName
                           << " hits were generated in the event."
                           << endl
                           << "Only " << _sizeLimit << " are saved in output collection."
                           << endl;

      G4cout << "Total of " << _currentSize << " "
           << SensitiveDetectorName
           << " hits were generated in the event."
           << G4endl
           << "Only " << _sizeLimit << " are saved in output collection."
           << G4endl;

    }

    if (verboseLevel>0) {
      G4int NbHits = _collection->size();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits
             << " hits in " << SensitiveDetectorName << ": " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i].print(G4cout);
    }

  }


  void Mu2eSensitiveDetector::beforeG4Event(StepPointMCCollection& outputHits, 
                                           PhysicsProcessInfo& processInfo,
                                           art::ProductID const& simID,
                                           art::Event const & event ){
    _collection  = &outputHits;
    _processInfo = &processInfo;
    _simID       = &simID;
    _event       = &event;

    return;

  }


} //namespace mu2e
