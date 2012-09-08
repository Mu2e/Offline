//
// Define a sensitive detector for CaloCrystal Detectors
//
// $Id: CaloCrystalSD.cc,v 1.20 2012/09/08 02:24:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/09/08 02:24:25 $
//
// Original author Ivan Logashenko
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/WorldG4.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  CaloCrystalSD::CaloCrystalSD(G4String name, SimpleConfig const & config ): 
    Mu2eSensitiveDetector(name,config)
  { }
  
  G4bool CaloCrystalSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    // Calculate energy deposition in the crystal
    G4double edep = aStep->GetTotalEnergyDeposit();
    if( edep<=0 ) return false;

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in " 
                             << SensitiveDetectorName
                             << ": "
                             << _currentSize << endl;
      }
      return false;
    }

    const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();

    // Get crystal ID
    G4int copyNo = touchableHandle->GetCopyNumber(0);

    G4String const& pname  = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    ProcessCode endCode(_processInfo->findAndCount(pname));

   // Originally the hit position was saved in local crystal frame.
    // Not it is saved in Mu2e frame, hence the following code is
    // commented out.
    // Calculate enerdy deposition position along the crystal
     G4AffineTransform const& toLocal = touchableHandle->GetHistory()->GetTopTransform();
     G4AffineTransform        toWorld = toLocal.Inverse();
     G4ThreeVector posWorld = aStep->GetPreStepPoint()->GetPosition();
     G4ThreeVector posLocal = toLocal.TransformPoint(posWorld);


    // Add the hit to the framework collection.
    // The point's coordinates are saved in the mu2e coordinate system.
    _collection->
      push_back(StepPointMC(art::Ptr<SimParticle>
                            ( *_simID,
                              aStep->GetTrack()->GetTrackID(),
                              _event->productGetter(*_simID) ),
                            copyNo,
                            edep,
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

} //namespace mu2e
