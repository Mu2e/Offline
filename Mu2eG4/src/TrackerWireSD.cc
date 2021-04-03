//
// Sensitive detector for TrackerWire
//

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/TrackerWireSD.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"

// G4 includes
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4ios.hh"

using namespace std;

namespace mu2e {

  G4ThreeVector TrackerWireSD::_mu2eDetCenter;

  TrackerWireSD::TrackerWireSD(G4String name, const SimpleConfig& config) :
                  Mu2eG4SensitiveDetector(name,config) { }

  TrackerWireSD::~TrackerWireSD(){ }

  G4bool TrackerWireSD::ProcessHits(G4Step* aStep, G4TouchableHistory*){

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of steps reached "
                             << SensitiveDetectorName
                             << ": "
                             << _currentSize << endl;
      }
      return false;
    }

    G4double edep  = aStep->GetTotalEnergyDeposit();
    G4double nidep = aStep->GetNonIonizingEnergyDeposit();
    G4double stepL  = aStep->GetStepLength();
    //G4double idep  = edep-nidep;

    if ( _debugList.inList() )  G4cout<<"edep "<<edep<<" nidep "<<nidep<<" step "<<stepL<<G4endl;
    // I am not sure why we get these cases but we do.  Skip them.
    if ( (edep == 0. /*|| idep == 0.*/)/*&& step == 0.*/ ) {
      if ( _debugList.inList() )  G4cout<<"Skipped"<<G4endl;
      return false;
    }

    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4int motherCopyNo = preStepPoint->GetTouchableHandle()->GetReplicaNumber(1);

    if ( _debugList.inList() )  {
      G4cout<<"Step vol name "<<aStep->GetTrack()->GetVolume()->GetName()<<G4endl;
      G4cout<<"Step vol copyNo "<<" mother copyNo "<<motherCopyNo<<G4endl;
    }

    // Which process caused this step to end?
    ProcessCode endCode(_processInfo->
                        findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));

    // Add the hit to the framework collection.
    // The point's coordinates are saved in the mu2e coordinate system.
    _collection->
      push_back(StepPointMC(_spHelper->particlePtr(aStep->GetTrack()),
                            motherCopyNo,
                            edep,
                            nidep,
                            0., // visible energy deposit; used in scintillators
                            preStepPoint->GetGlobalTime(),
                            preStepPoint->GetProperTime(),
                            preStepPoint->GetPosition() - _mu2eDetCenter,
                            aStep->GetPostStepPoint()->GetPosition() - _mu2eDetCenter,
                            preStepPoint->GetMomentum(),
                            aStep->GetPostStepPoint()->GetMomentum(),
                            stepL,
                            endCode
                            ));

    return true;
  }

} //namespace mu2e
