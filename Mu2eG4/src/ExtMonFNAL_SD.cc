#include <cstdio>

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4ios.hh"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/ExtMonFNAL_SD.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"

#define AGDEBUG(stuff) std::cerr<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;

using namespace std;

namespace mu2e {
  const int EMFrameToSensorDepth = 1;   // The sensitive volume is the direct daughter of the volume that defines the EMF frame

  ExtMonFNAL_SD::ExtMonFNAL_SD(G4String name, const SimpleConfig& config) :
    G4VSensitiveDetector(name),
    _collection(0),
    _processInfo(0),
    _mu2eOrigin(GeomHandle<WorldG4>()->mu2eOriginInWorld()),
    _sizeLimit(config.getInt("g4.stepsSizeLimit",0)),
    _currentSize(0),
    _simID(0),
    _event(0)
  {
  }

  void ExtMonFNAL_SD::Initialize(G4HCofThisEvent* HCE){
    _currentSize=0;
  }


  G4bool ExtMonFNAL_SD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    // Calculate energy deposition
    G4double edep = aStep->GetTotalEnergyDeposit();

    AGDEBUG("AG: edep = "<<edep);


    if( edep<=0 ) return false;

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in ExtMonFNAL_SD: "
                             << _currentSize << endl;
      }
      return false;
    }

    //G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    //G4int eventId = event->GetEventID();
    //G4int trackId = aStep->GetTrack()->GetTrackID();

    // These are in G4 world coordinates
    G4ThreeVector prePosG4 = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector preMomG4 = aStep->GetPreStepPoint()->GetMomentum();

    AGDEBUG("AG: got prePosG4     = "<<prePosG4<<    ", preMomG4     = "<<preMomG4);

    // Transform the global coordinates into the ExtMon frame
    const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();
    const G4AffineTransform& toEMFrame = touchableHandle->GetHistory()->GetTransform(EMFrameToSensorDepth);
    G4ThreeVector prePosExtMon = toEMFrame.TransformPoint(prePosG4);
    G4ThreeVector preMomExtMon = toEMFrame.TransformAxis(preMomG4);

    // FIXME: add a check comparing transormation results using G4 and using the ExtMonFNAL geometry class.


    
    AGDEBUG("AG: got prePosExtMon = "<<prePosExtMon<<", preMomExtMon = "<<preMomExtMon);

    //: const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();
    //: static G4ThreeVector detectorOrigin = GetTrackerOrigin(touchableHandle);
    //: // Position at start of step point, in world system and in
    //: // a system in which the center of the tracking detector is the origin.
    //: G4ThreeVector prePosWorld   = aStep->GetPreStepPoint()->GetPosition();
    //: G4ThreeVector prePosTracker = prePosWorld - detectorOrigin;
    //: G4ThreeVector preMomWorld = aStep->GetPreStepPoint()->GetMomentum();

    //: // from http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch04.html#sect.Geom.Touch
    //: 
    //: G4ThreeVector worldPosition = preStepPoint->GetPosition();
    //: G4ThreeVector localPosition = theTouchable->GetHistory()->
    //:   GetTopTransform().TransformPoint(worldPosition);



    // Which process caused this step to end?
    G4String const& pname  = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    AGDEBUG("pname = "<<pname);
    ProcessCode endCode(_processInfo->findAndCount(pname));

    unsigned pixelSensorID = 0; // AG: FIXME: map it from the volume pointer ...;

    // Add the hit to the framework collection.
    // The point's coordinates are saved in the mu2e coordinate system.
    _collection->
      push_back(StepPointMC(art::Ptr<SimParticle>( *_simID, aStep->GetTrack()->GetTrackID(), _event->productGetter(*_simID) ),
                            pixelSensorID,
                            edep,
                            aStep->GetNonIonizingEnergyDeposit(),
                            aStep->GetPreStepPoint()->GetGlobalTime(),
                            aStep->GetPreStepPoint()->GetProperTime(),
			    prePosExtMon,
                            preMomExtMon,
                            aStep->GetStepLength(),
                            endCode
                            ));

    return true;

  }


  void ExtMonFNAL_SD::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      mf::LogWarning("G4") << "Total of " << _currentSize
                           << " calorimeter hits were generated in the event."
                           << endl
                           << "Only " << _sizeLimit << " are saved in output collection."
                           << endl;
    }

    if (verboseLevel>0) {
      G4int NbHits = _collection->size();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits
             << " hits in the calorimeter: " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i].print(G4cout);
    }

  }


  void ExtMonFNAL_SD::beforeG4Event(StepPointMCCollection& outputHits,
                                    PhysicsProcessInfo& processInfo,
                                    art::ProductID const& simID,
                                    art::Event const & event ){
    _collection    = &outputHits;
    _processInfo   = &processInfo;
    _simID         = &simID;
    _event         = &event;

    return;
  } // end of beforeG4Event


} //namespace mu2e
