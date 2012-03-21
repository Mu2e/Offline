//
// Define a sensitive detector for TTrackerDeviceSupport
//
// $Id: TTrackerDeviceSupportSD.cc,v 1.2 2012/03/21 15:52:09 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/21 15:52:09 $
//
// Original author KLG
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/TTrackerDeviceSupportSD.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "TTrackerGeom/inc/TTracker.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  TTrackerDeviceSupportSD::TTrackerDeviceSupportSD(G4String name, const SimpleConfig& config):
    G4VSensitiveDetector(name),
    _collection(0),
    _processInfo(0),
    _mu2eOrigin(GeomHandle<WorldG4>()->mu2eOriginInWorld()),
    _TrackerVersion(0),
    _debugList(0),
    _sizeLimit(config.getInt("g4.stepsSizeLimit",0)),
    _currentSize(0),
    _simID(0),
    _event(0)
  {
    SetVerboseLevel(config.getInt("ttracker.verbosityLevel",0));

    // Get list of events for which to make debug printout.
    string key("g4.ttrackerDeviceSupportSDEventList");
    if ( config.hasName(key) ){
      vector<int> list;
      config.getVectorInt(key,list);
      _debugList.add(list);
    }

   art::ServiceHandle<GeometryService> geom;

    if ( !geom->hasElement<TTracker>() ) {
      throw cet::exception("GEOM")
        << "Expected T Tracker but did not find it.\n";
    } 
    else {
      _TrackerVersion = config.getInt("TTrackerVersion",3);
      if ( _TrackerVersion != 3) {
        throw cet::exception("TTrackerDeviceSupportSD")
          << "Expected TTrackerVersion of 3 but found " << _TrackerVersion <<endl;
      }
    }
  }


  TTrackerDeviceSupportSD::~TTrackerDeviceSupportSD(){ }

  void TTrackerDeviceSupportSD::Initialize(G4HCofThisEvent* HCE){

    _currentSize=0;

  }


  G4bool TTrackerDeviceSupportSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in TTrackerDeviceSupportSD: "
                             << _currentSize << endl;
      }
      return false;
    }

    // Which process caused this step to end?
    G4String const& pname  = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    ProcessCode endCode(_processInfo->findAndCount(pname));

    G4int sdcn = 0;

    G4TouchableHandle const & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();

    G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();


    // we take here the numbering from TTrackerDeviceEnvelope
    // (not TTrackerDeviceSupport which is always 0)

    G4int en = event->GetEventID();
    G4int ti = aStep->GetTrack()->GetTrackID();
    G4int cn = touchableHandle->GetCopyNumber(1);

    //    G4TouchableHistory* theTouchable =
    //  (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );

    if (verboseLevel > 1) {

      cout << "TTrackerDeviceSupportSD::" << __func__ << " Debugging history depth " <<
        setw(4) << touchableHandle->GetHistoryDepth() << endl;

      cout << "TTrackerDeviceSupportSD::" << __func__ << " Debugging copy n 0 1 2 3 4 " <<
        setw(4) << touchableHandle->GetCopyNumber(0) <<
        setw(4) << touchableHandle->GetCopyNumber(1) <<
        setw(4) << touchableHandle->GetCopyNumber(2) <<
        setw(4) << touchableHandle->GetCopyNumber(3) <<
        setw(4) << touchableHandle->GetCopyNumber(4) << endl;

      cout << "TTrackerDeviceSupportSD::" << __func__ << " Debugging PV Name Mother Name " <<
        touchableHandle->GetVolume(0)->GetName() << " " <<
        touchableHandle->GetVolume(1)->GetName() << " " <<
        touchableHandle->GetVolume(2)->GetName() << " " <<
        touchableHandle->GetVolume(3)->GetName() << " " <<
        touchableHandle->GetVolume(4)->GetName() << endl;

      cout << "TTrackerDeviceSupportSD::" << __func__ << " Debugging hit info event track copyn replican: " <<
        setw(4) << en << " " <<
        setw(4) << ti << " " <<
        setw(4) << cn << endl;

      cout << "TTrackerDeviceSupportSD::" << __func__ << " Debugging _TrackerVersion: " << 
        _TrackerVersion << endl;

    }
    // device support is placed in a device envelope; device name should be unique though
    // therefore _TrackerVersion is probably not needed
    //

    if ( _TrackerVersion == 3) {
      sdcn = cn;
    }

    // Add the hit to the framework collection.
    // The point's coordinates are saved in the mu2e coordinate system.
    _collection->
      push_back(StepPointMC(art::Ptr<SimParticle>( *_simID, aStep->GetTrack()->GetTrackID(), _event->productGetter(*_simID) ),
                            sdcn,
                            aStep->GetTotalEnergyDeposit(),
                            aStep->GetNonIonizingEnergyDeposit(),
                            aStep->GetPreStepPoint()->GetGlobalTime(),
                            aStep->GetPreStepPoint()->GetProperTime(),
                            aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin,
                            aStep->GetPreStepPoint()->GetMomentum(),
                            aStep->GetStepLength(),
                            endCode
                            ));

    if (verboseLevel >0) {
      cout << "TTrackerDeviceSupportSD::" << __func__ << " Event " << setw(4) <<
        G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() <<
        " TTrackerDeviceSupport " << 
        touchableHandle->GetVolume(1)->GetName() << " " <<
        setw(4) << touchableHandle->GetVolume(1)->GetCopyNo() <<
        " hit at: " << aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin << endl;
    }

    // Debugging tests (empty for now)
    if ( !_debugList.inList() ) return true;
    return true;

  }


  void TTrackerDeviceSupportSD::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      mf::LogWarning("G4") << "Total of " << _currentSize
                           << " TTrackerDeviceSupport detector hits were generated in the event."
                           << endl
                           << "Only " << _sizeLimit << " are saved in output collection."
                           << endl;
      cout << "Total of " << _currentSize
           << " TTrackerDeviceSupport detector hits were generated in the event."
           << endl
           << "Only " << _sizeLimit << " are saved in output collection."
           << endl;
    }

    if (verboseLevel>0) {
      G4int NbHits = _collection->size();
      G4cout << "\n-------->Hit Collection: in this event there are " << NbHits
             << " hits in the TTrackerDeviceSupport: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i].print(G4cout);
    }

  }


  void TTrackerDeviceSupportSD::beforeG4Event(StepPointMCCollection& outputHits, 
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
