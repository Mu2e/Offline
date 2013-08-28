//
// Define a sensitive detector for extinction monitor tof
//
// $Id: ExtMonUCITofSD.cc,v 1.5 2013/08/28 05:58:17 gandr Exp $
// $Author: gandr $
// $Date: 2013/08/28 05:58:17 $
//

#include <cstdio>
#include <fstream>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4ios.hh"

// Mu2e includes
#include "Mu2eG4/inc/ExtMonUCITofSD.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"

using namespace std;

namespace mu2e {

  ExtMonUCITofSD::ExtMonUCITofSD(G4String name, const SimpleConfig& config):
    Mu2eSensitiveDetector(name,config)
  { }

  G4bool ExtMonUCITofSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in ExtMonUCITofSD: "
                              << _currentSize << endl;
      }
      return false;
    }

    //G4Event const* event = G4RunManager::GetRunManager()->GetCurrentEvent();
    //const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();
    //G4int eventId = event->GetEventID();
    //G4int trackId = aStep->GetTrack()->GetTrackID();
    //G4int copyNo = touchableHandle->GetCopyNumber();

    // Get tof ID
    //G4int copyNo = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2);

    // Which process caused this step to end?
    ProcessCode endCode(_processInfo->
                        findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));

    // Add the hit to the framework collection.
    // The point's coordinates are saved in the mu2e coordinate system.
    _collection->
      push_back(StepPointMC(_spHelper->particlePtr(aStep->GetTrack()),
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

      int static const verbosityLevel = 0;
      if (verbosityLevel >0) {
       cout << __func__ << " Event " << 
         G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID() <<
         " ExtMonTof " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << " " <<
         aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetCopyNo() <<
         " hit at: " << aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin << endl;
       }

    return true;

  }

} //namespace mu2e
