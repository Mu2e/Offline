//
// Sensitive detector for TrackerPlaneSupport
//
// Original author KLG
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/TrackerPlaneSupportSD.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"

// G4 includes
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4ios.hh"

using namespace std;

namespace mu2e {

  TrackerPlaneSupportSD::TrackerPlaneSupportSD(G4String name, SimpleConfig const & config ):
    Mu2eG4SensitiveDetector(name,config)
  {

    SetVerboseLevel(config.getInt("tracker.verbosityLevel",0));

   art::ServiceHandle<GeometryService> geom;

    if ( !geom->hasElement<Tracker>() ) {
      throw cet::exception("GEOM")
        << "Expected Tracker but did not find it.\n";
    }
    else {
      _TrackerVersion = config.getInt("TrackerVersion",3);
      if ( _TrackerVersion < 3) {
        throw cet::exception("TrackerPlaneSupportSD")
          << "Expected TrackerVersion >= 3 but found " << _TrackerVersion <<endl;
      }
    }
  }

  G4bool TrackerPlaneSupportSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    _currentSize += 1;

    if ( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of steps reached in "
                             << SensitiveDetectorName
                             << ": "
                             << _currentSize << endl;
      }
      return false;
    }

    // Which process caused this step to end?
    ProcessCode endCode(_processInfo->
                        findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));

    G4int sdcn = 0;

    G4TouchableHandle const & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();

    G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();


    // we take here the numbering from TrackerPlaneEnvelope
    // (not TrackerPlaneSupport which is always 0)

    G4int en = event->GetEventID();
    G4int ti = aStep->GetTrack()->GetTrackID();
    G4int cn = touchableHandle->GetCopyNumber(1);

    //    G4TouchableHistory* theTouchable =
    //  (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );

    if (verboseLevel > 1) {

      cout << "TrackerPlaneSupportSD::" << __func__ << " Debugging history depth " <<
        setw(4) << touchableHandle->GetHistoryDepth() << endl;

      cout << "TrackerPlaneSupportSD::" << __func__ << " Debugging copy n 0 1 2 3 4 " <<
        setw(4) << touchableHandle->GetCopyNumber(0) <<
        setw(4) << touchableHandle->GetCopyNumber(1) <<
        setw(4) << touchableHandle->GetCopyNumber(2) <<
        setw(4) << touchableHandle->GetCopyNumber(3) <<
        setw(4) << touchableHandle->GetCopyNumber(4) << endl;

      cout << "TrackerPlaneSupportSD::" << __func__ << " Debugging PV Name Mother Name " <<
        touchableHandle->GetVolume(0)->GetName() << " " <<
        touchableHandle->GetVolume(1)->GetName() << " " <<
        touchableHandle->GetVolume(2)->GetName() << " " <<
        touchableHandle->GetVolume(3)->GetName() << " " <<
        touchableHandle->GetVolume(4)->GetName() << endl;

      cout << "TrackerPlaneSupportSD::" << __func__
           << " Debugging hit info event track copyn replican: " <<
        setw(4) << en << " " <<
        setw(4) << ti << " " <<
        setw(4) << cn << endl;

      cout << "TrackerPlaneSupportSD::" << __func__ << " Debugging _TrackerVersion: " <<
        _TrackerVersion << endl;

    }
    // plane support is placed in a plane envelope; plane name should be unique though
    // therefore _TrackerVersion is probably not needed
    //

    if ( _TrackerVersion == 3) {
      sdcn = cn;
    }

    // Add the hit to the framework collection.
    // The point's coordinates are saved in the mu2e coordinate system.
    _collection->
      push_back(StepPointMC(_spHelper->particlePtr(aStep->GetTrack()),
                            sdcn,
                            aStep->GetTotalEnergyDeposit(),
                            aStep->GetNonIonizingEnergyDeposit(),
                            0., // visible energy deposit; used in scintillators
                            aStep->GetPreStepPoint()->GetGlobalTime(),
                            aStep->GetPreStepPoint()->GetProperTime(),
                            aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin,
                            aStep->GetPostStepPoint()->GetPosition() - _mu2eOrigin,
                            aStep->GetPreStepPoint()->GetMomentum(),
                            aStep->GetPostStepPoint()->GetMomentum(),
                            aStep->GetStepLength(),
                            endCode
                            ));

    if (verboseLevel >0) {
      cout << "TrackerPlaneSupportSD::" << __func__ << " Event " << setw(4) <<
        G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() <<
        " TrackerPlaneSupport " <<
        touchableHandle->GetVolume(1)->GetName() << " " <<
        setw(4) << touchableHandle->GetVolume(1)->GetCopyNo() <<
        " hit at: " << aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin << endl;
    }

    // Debugging tests (empty for now)
    if ( !_debugList.inList() ) return true;
    return true;

  }

} //namespace mu2e
