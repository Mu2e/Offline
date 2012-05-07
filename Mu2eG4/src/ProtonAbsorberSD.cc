//
// Define a sensitive detector for proton absorber
//
//  $Id: ProtonAbsorberSD.cc,v 1.1 2012/05/07 23:35:58 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/05/07 23:35:58 $
//
// Original author MyeongJae Lee
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/ProtonAbsorberSD.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  ProtonAbsorberSD::ProtonAbsorberSD(G4String name, const SimpleConfig& config):
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
    string key("g4.virtualSDEventList");
    if ( config.hasName(key) ){
      vector<int> list;
      config.getVectorInt(key,list);
      _debugList.add(list);
    }

  }


  ProtonAbsorberSD::~ProtonAbsorberSD(){ }

  void ProtonAbsorberSD::Initialize(G4HCofThisEvent* HCE){

    _currentSize=0;

  }


  G4bool ProtonAbsorberSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in ProtonAbsorberSD: "
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
      push_back(StepPointMC(art::Ptr<SimParticle>( *_simID, aStep->GetTrack()->GetTrackID(), _event->productGetter(*_simID) ),
                            aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetCopyNo(),
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


  void ProtonAbsorberSD::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      mf::LogWarning("G4") << "Total of " << _currentSize
                           << " proton absorber hits were generated in the event."
                           << endl
                           << "Only " << _sizeLimit << " are saved in output collection."
                           << endl;
      cout << "Total of " << _currentSize
           << " proton absorber hits were generated in the event."
           << endl
           << "Only " << _sizeLimit << " are saved in output collection."
           << endl;
    }

    if (verboseLevel>0) {
      G4int NbHits = _collection->size();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits
             << " hits in the proton absorber: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i].print(G4cout);
    }

  }


  void ProtonAbsorberSD::beforeG4Event(StepPointMCCollection& outputHits, 
                                       PhysicsProcessInfo& processInfo,
                                       art::ProductID const& simID,
                                       art::Event const & event ){
    _collection    = &outputHits;
    _processInfo   = &processInfo;
    _simID         = &simID;
    _event = &event;

    return;
  } // end of beforeG4Event


} //namespace mu2e
