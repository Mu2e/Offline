//
// DefinesGeant4  sensitive detector for a typicaly numbered volume using Mu2e reference frame
//
// Original author KLG
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4SensitiveDetector.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"

// G4 includes
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4ios.hh"
#include "Geant4/G4Threading.hh"

using namespace std;

namespace mu2e {

  Mu2eG4SensitiveDetector::Mu2eG4SensitiveDetector(G4String const name, const SimpleConfig& config):
    G4VSensitiveDetector(name),
    _collection(0),
    _processInfo(0),
    _mu2eOrigin(GeomHandle<WorldG4>()->mu2eOriginInWorld()),
    _debugList(0),
    _sizeLimit(config.getInt("g4.stepsSizeLimit",0)),
    _currentSize(0),
    _spHelper()
  {

   // Get list of events for which to make debug printout.
    // we generate the key string based on the detector name
    // consult Mu2eG4/inc/SensitiveDetectorName.hh for the names
    std::ostringstream sdKeyName;
    sdKeyName<<"g4."<< SensitiveDetectorName << "SDEventList";
    // G4cout << __func__ << " sdKeyName: " << sdKeyName.str() << G4endl;
    // G4cout << __func__ << " sd name: " << name << G4endl;

    string key(sdKeyName.str());
    if ( config.hasName(key) ){
      vector<int> list;
      config.getVectorInt(key,list);
      _debugList.add(list);
    }
  }

  void Mu2eG4SensitiveDetector::Initialize(G4HCofThisEvent* HCE){
      _currentSize=0;
  }


  G4bool Mu2eG4SensitiveDetector::ProcessHits(G4Step* aStep,G4TouchableHistory*){

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

// this little section of code containing 'if ( _debugList.inList() )'
// was occasionally causing seg faults in MT mode
// this issue seems to have been fixed
// but I have taken it out to reduce output 05/18

/*    if ( _debugList.inList() )  {
            G4cout<<"edep "<<aStep->GetTotalEnergyDeposit()
                  <<" nidep "<<aStep->GetNonIonizingEnergyDeposit()
                  <<" step "<<aStep->GetStepLength()<<G4endl;
            G4cout<<"Step vol name "<<aStep->GetTrack()->GetVolume()->GetName()<<G4endl;
    }
*/

    // Which process caused this step to end?
    ProcessCode endCode(_processInfo->
                findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));

      // Add the hit to the framework collection.
      // The point's coordinates are saved in the mu2e coordinate system.
    _collection->
      push_back(StepPointMC(_spHelper->particlePtr(aStep->GetTrack()),
                            aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(),
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
      return true;

  }//ProcessHits


  void Mu2eG4SensitiveDetector::EndOfEvent(G4HCofThisEvent*){

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
      for (G4int i=0;i<NbHits;i++) (*_collection)[i].print(G4cout, true, false);
    }

  }//EndOfEvent


  void Mu2eG4SensitiveDetector::beforeG4Event(StepPointMCCollection& outputHits,
                                            PhysicsProcessInfo& processInfo,
                                            const SimParticleHelper& spHelper){
    _collection  = &outputHits;
    _processInfo = &processInfo;
    _spHelper    = &spHelper;

    return;

  }//beforeG4Event

} //namespace mu2e
